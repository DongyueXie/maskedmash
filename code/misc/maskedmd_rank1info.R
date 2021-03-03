#'@title masked deconvolution: estimate prior covariance matrices using multivariate deconvolution
#'@description This function estimates prior weights and variances usnig masked z scores. If strong
#'is provided, then U.canon should be set to NULL and usepointmass set to FALSE;
#'@param strong index of strong effect samples.
#'
masked.md.rank1info = function(data,strong=NULL,thresh=NULL,
                     pi=NULL,U.canon=NULL,U.data,
                     usepointmass=FALSE,
                     pi_thresh = 1e-8,
                     max_iter=1000,
                     tol=1e-3,
                     verbose=TRUE,
                     printevery = 1,
                     mc.cores = 4){

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

  lambda = c()
  U = matrix(nrow=R,ncol = length(U.data))

  t.lik = 0
  t.obj = 0
  t.postw = 0
  t.pi = 0
  t.U = 0


  for(iter in 1:max_iter){



    t1 = Sys.time()
    # calculate marginal likelihood
    # a list of length N, each a J*K log likelihood matrix
    lik.list = lapply(Z.comb,function(x){
      temp = calc_relative_lik_matrix(mash_set_data(x$z.comb),Ulist,algorithm.version = 'Rcpp')
      temp
    })

    t2 = Sys.time()
    t.lik = t.lik + t2 - t1

    # evaluate objective function
    # compute_vloglik_from_matrix_and_pi returns a J*1 loglik vector(after summing over k)
    loglik.obj[iter+1] = sum(log(unlist(lapply(1:N,function(i){
      lm_i = lik.list[[i]]
      det_i = Z.comb[[i]]$z.det.jacob
      sum(det_i*exp(compute_vloglik_from_matrix_and_pi(pi,lm_i, matrix(1,nrow=length(det_i),ncol=R))))
    }))))

    t3 = Sys.time()
    t.obj = t.obj + t3-t2

    ## check convergence
    if(abs(loglik.obj[iter+1]-loglik.obj[iter])<tol){
      break
    }

    # update gamma_ijk
    post_weights = calc_post_weights(Z.comb,lik.list,pi)

    t4 = Sys.time()
    t.postw = t.postw + t4-t3

    ## update pi
    pi = lapply(post_weights,colSums)
    pi = Reduce('+',pi)/N

    t5 = Sys.time()
    t.pi = t.pi + t5-t4

    # update U.data
    ## eigenvectors of sum_ij(gamma_ijk S_ij)
    ## let gamma.tilde_ij = sum_l gamma_ijkl w_l/(1+w_l)

    # update.Uk = mclapply((usepointmass+length(U.canon)+1):K,function(k){
    #   S.tilde_k = 0
    #   for(i in 1:N){
    #     Ji = nrow(Z.comb[[i]]$z.comb)
    #     for(j in 1:Ji){
    #       S.tilde_k = S.tilde_k + post_weights[[i]][j,k]*tcrossprod(Z.comb[[i]]$z.comb[j,])
    #     }
    #   }
    #   temp = eigen(S.tilde_k/(N*pi[k]))
    #
    #   u_k = temp$vectors[,1]
    #   U[,(k-usepointmass-length(U.canon))] = u_k
    #   lambda_k = max(temp$values[1]-1, 0)
    #   lambda[(k-usepointmass-length(U.canon))] = lambda_k
    #   Ulist[[k]] = lambda_k*tcrossprod(u_k)
    # })

    for(k in (usepointmass+length(U.canon)+1):K){

      # S.tilde_k = 0
      # for(i in 1:N){
      #   Ji = nrow(Z.comb[[i]]$z.comb)
      #   for(j in 1:Ji){
      #     S.tilde_k = S.tilde_k + post_weights[[i]][j,k]*tcrossprod(Z.comb[[i]]$z.comb[j,])
      #   }
      # }

      S.tilde_k = lapply(1:N,function(i){
        crossprod(sqrt(post_weights[[i]][,k])*Z.comb[[i]]$z.comb)
      })

      S.tilde_k = Reduce("+",S.tilde_k)/(N*pi[k])

      temp = eigen(S.tilde_k)

      u_k = temp$vectors[,1]
      U[,(k-usepointmass-length(U.canon))] = u_k
      lambda_k = max(temp$values[1]-1, 0)
      lambda[(k-usepointmass-length(U.canon))] = lambda_k
      Ulist[[k]] = lambda_k*tcrossprod(u_k)

      #Ulist[[k]] = temp$vectors%*%diag(pmax(temp$values-1,0))%*%t(temp$vectors)
      #Ulist[[k]] = temp$vectors%*%diag(pmax(temp$values,0))%*%t(temp$vectors)

    }

    t6 = Sys.time()
    t.U = t.U + t6-t5

    if(verbose){
      if(iter%%printevery==0){
        print(sprintf("Done %1.0f iterations, loglik at %.2f",iter,loglik.obj[iter+1]))
      }
    }

  }

  U.est = Ulist[(usepointmass+length(U.canon)+1):K]
  U.est.adj = U.est

  # calculate second derivatives
  Hessians = list()
  for(k in 1:length(U.data)){
    Hessian = 0
    for(i in 1:N){
      Ji = nrow(Z.comb[[i]]$z.comb)
      for(j in 1:Ji){
        Hessian = Hessian - post_weights[[i]][j,k]*(tcrossprod(Z.comb[[i]]$z.comb[j,])/(1+lambda[k]))
      }
    }
    Hessians[[k]] = Hessian
    U.est.adj[[k]] = tcrossprod(U[,k]) + solve(-Hessian)
  }


  #################### calculate information matrixm#################

  #browser()
  # Inv.info = list()
  # Hessians = list()
  # for(k in 1:length(U.data)){
  #
  #   Hessian = 0
  #   Ell = 0
  #   El = 0
  #
  #   for(i in 1:N){
  #     Ji = nrow(Z.comb[[i]]$z.comb)
  #     H_i = matrix(0,nrow=R+2,ncol=R+2)
  #     l_ij = c()
  #     Ell_i = 0
  #     #El_i = 0
  #     for(j in 1:Ji){
  #       H_i[1,1] = H_i[1,1] - post_weights[[i]][j,k]/pi[k]^2
  #       H_i[2,2] = H_i[2,2] + post_weights[[i]][j,k]*(1/(2*(1+lambda[k])^2)-
  #                                                       t(U[,k])%*%tcrossprod(Z.comb[[i]]$z.comb[j,])%*%U[,k]/(1+lambda[k])^3)
  #       #H_i[3:(R+2),3:(R+2)] = H_i[3:(R+2),3:(R+2)] + post_weights[[i]][j,k]*(tcrossprod(Z.comb[[i]]$z.comb[j,])*lambda[k]/(1+lambda[k]))
  #       H_i[3:(R+2),3:(R+2)] = H_i[3:(R+2),3:(R+2)] - post_weights[[i]][j,k]*(tcrossprod(Z.comb[[i]]$z.comb[j,])/(1+lambda[k]))
  #       H_i[3:(R+2),2] = H_i[3:(R+2),2] + post_weights[[i]][j,k]*(1/(1+lambda[k])^2*tcrossprod(Z.comb[[i]]$z.comb[j,])%*%U[,k])
  #       H_i[2,3:(R+2)] = H_i[3:(R+2),2]
  #
  #       l_ij[1] = 1/pi[k]
  #       l_ij[2] = -1/(2*(1+lambda[k])) + t(U[,k])%*%tcrossprod(Z.comb[[i]]$z.comb[j,])%*%U[,k]/(2*(1+lambda[k])^2)
  #       #l_ij[3:(R+2)] = lambda[k]/(1+lambda[k])*tcrossprod(Z.comb[[i]]$z.comb[j,])%*%U[,k]
  #       l_ij[3:(R+2)] = -1/(1+lambda[k])*tcrossprod(Z.comb[[i]]$z.comb[j,])%*%U[,k]
  #       Ell_i = Ell_i + post_weights[[i]][j,k]*tcrossprod(l_ij)
  #
  #       #El_i = El_i +
  #       El = El+ post_weights[[i]][j,k]*l_ij
  #
  #
  #       #Ell = Ell + post_weights[[i]][j,k]*tcrossprod(lambda[k]/(1+lambda[k])*(tcrossprod(Z.comb[[i]]$z.comb[j,])%*%U[,k]))
  #       #El = El + post_weights[[i]][j,k]*(lambda[k]/(1+lambda[k])*(tcrossprod(Z.comb[[i]]$z.comb[j,])%*%U[,k]))
  #     }
  #     Hessian = Hessian+H_i
  #     Hessians[[k]] = Hessian
  #     Ell = Ell+Ell_i
  #
  #   }
  #
  #   #browser()
  #   I_k = -Hessian - Ell + tcrossprod(El)
  #   Inv.info[[k]] = solve(I_k)



  #   }


  ## now we have estimated prior weights and covariance
  ## calculate posterior summaries

  #result = calc_post_summary(Z.comb,xUlist,pi,post_weights)

  t = list(t.lik = t.lik,t.obj=t.obj,t.postw=t.postw,t.pi=t.pi,t.U=t.U)


  result = list()
  result$U.est = U.est
  result$pi = pi
  result$loglik = loglik.obj
  result$U = U
  result$U.est.adj = U.est.adj
  result$Hessians = Hessians
  result$lambda = lambda
  result$post_weights = post_weights
  result$t = t

  result
}
