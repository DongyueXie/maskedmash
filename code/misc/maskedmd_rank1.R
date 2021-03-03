
#'@title masked deconvolution: estimate prior covariance matrices(rank 1) using multivariate deconvolution
#'@description This function estimates prior weights and variances usnig masked z scores. If strong
#'is provided, then U.canon should be set to NULL and usepointmass set to FALSE;
#'@param strong index of strong effect samples.
#'
masked.md.rank1 = function(data,strong=NULL,thresh=NULL,
                     pi=NULL,U.canon=NULL,U.data,
                     usepointmass=FALSE,
                     pi_thresh = 1e-8,
                     max_iter=1000,
                     tol=1e-3,
                     verbose=TRUE,
                     printevery = 1){

  Z = data$Bhat
  if(is.null(strong)){
    strong = 1:nrow(Z)
  }
  Z = Z[strong,]
  N = nrow(Z)
  R = ncol(Z)
  I_R = diag(R)

  if(is.null(thresh)){
    # mask all a-scores
    thresh = qnorm(0.75)
  }

  Z.comb = set_Z.comb(Z,thresh)
  Ulist = c(U.canon,U.data)
  if(usepointmass){
    Ulist = c(list(null=matrix(0,nrow = R,ncol = R)),Ulist)
  }else{
    Ulist = Ulist
  }
  K = length(Ulist)

  if(is.null(pi)){pi = rep(1/K,K)}
  loglik.obj = -Inf

  theta = c()
  lambda = c()

  for(iter in 1:max_iter){



    # calculate marginal likelihood
    # a list of length N, each a J*K log likelihood matrix
    lik.list = lapply(Z.comb,function(x){
      temp = calc_relative_lik_matrix(mash_set_data(x$z.comb),Ulist)
      temp
    })

    # evaluate objective function
    # compute_vloglik_from_matrix_and_pi returns a J*1 loglik vector(after summing over k)
    loglik.obj[iter+1] = sum(log(unlist(lapply(1:N,function(i){
      lm_i = lik.list[[i]]
      det_i = Z.comb[[i]]$z.det.jacob
      sum(det_i*exp(compute_vloglik_from_matrix_and_pi(pi,lm_i, matrix(1,nrow=length(det_i),ncol=R))))
    }))))


    ## check convergence
    if(abs(loglik.obj[iter+1]-loglik.obj[iter])<tol){
      break
    }

    # update gamma_ijk
    post_weights = calc_post_weights(Z.comb,lik.list,pi)

    ## update pi
    pi = lapply(post_weights,colSums)
    pi = Reduce('+',pi)/N

    # update U.data
    ## eigenvectors of sum_ij(gamma_ijk S_ij)
    ## let gamma.tilde_ij = sum_l gamma_ijkl w_l/(1+w_l)

    for(k in (usepointmass+length(U.canon)+1):K){

      S.tilde_k = 0
      for(i in 1:N){
        Ji = nrow(Z.comb[[i]]$z.comb)
        for(j in 1:Ji){
          S.tilde_k = S.tilde_k + post_weights[[i]][j,k]*tcrossprod(Z.comb[[i]]$z.comb[j,])
        }
      }
      temp = eigen(S.tilde_k/(N*pi[k]))

      u_k = temp$vectors[,1]
      theta_k = max(sum(temp$values[-1])/(R-1) - 1, 0)
      theta[(k-usepointmass-length(U.canon))] = theta_k
      lambda_k = max(temp$values[1] - sum(temp$values[-1])/(R-1), 0)
      lambda[(k-usepointmass-length(U.canon))] = lambda_k
      Ulist[[k]] = theta_k*diag(R) + lambda_k*tcrossprod(u_k)
      #Ulist[[k]] = temp$vectors%*%diag(pmax(temp$values,0))%*%t(temp$vectors)

    }

    if(verbose){
      if(iter%%printevery==0){
        print(sprintf("Done %1.0f iterations, loglikelihood at %.2f",iter,loglik.obj[iter+1]))
      }
    }

  }

  ## now we have estimated prior weights and covariance
  ## calculate posterior summaries

  #result = calc_post_summary(Z.comb,xUlist,pi,post_weights)

  U.est = Ulist[(usepointmass+length(U.canon)+1):K]

  result = list()
  result$U.est = U.est
  result$pi = pi
  result$loglik = loglik.obj
  result$theta = theta
  result$lambda = lambda

  result
}
