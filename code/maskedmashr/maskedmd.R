#'@title masked deconvolution: estimate prior covariance matrices using multivariate deconvolution
#'@description This function estimates prior weights and variances usnig masked z scores. If strong
#'is provided, then U.canon should be set to NULL and usepointmass set to FALSE;
#'@param data object from mash_set_data
#'@param thresh NULL or a number, |z-score| larger than thresh will be masked.
#'@param pi prior weights; either fixed or as init value
#'@param U.canon a list of canonical prior cov matrices; fixed.
#'@param U.data a list of data-driven prior cov marrices; either fixed or as initialization.
#'@param algorithm.version either Rcpp or R
#'@param nu prior df of U. If NULL, nu = R+1
masked.md = function(data,strong=NULL,thresh=NULL,
                     pi=NULL,U.canon=NULL,U.data,
                     usepointmass=FALSE,
                     algorithm.version = 'Rcpp',
                     adjust = 'lb',
                     nu = NULL,
                     pi_thresh = 1e-8,
                     max_iter=1000,
                     tol=1e-3,
                     verbose=TRUE,
                     printevery = 50){

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

  Z.comb = mask.Z(Z,thresh)
  Ulist = c(U.canon,U.data)
  if(usepointmass){
    Ulist = c(list(null=matrix(0,nrow = R,ncol = R)),Ulist)
  }else{
    Ulist = Ulist
  }
  K = length(Ulist)

  if(is.null(pi)){pi = rep(1/K,K)}
  loglik.obj = -Inf

  for(iter in 1:max_iter){



    # calculate marginal likelihood
    # a list of length N, each a J*K log likelihood matrix
    lik.list = lapply(Z.comb,function(x){
      temp = calc_relative_lik_matrix(mash_set_data(x$z.comb),Ulist,algorithm.version = algorithm.version)
      temp
    })

    # evaluate objective function
    # compute_vloglik_from_matrix_and_pi returns a J*1 loglik vector(after summing over k)
    loglik.obj[iter+1] = calc_obj(lik.list,Z.comb,pi,N,R)


    #browser()
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

      S.tilde.k = lapply(1:N,function(i){
        crossprod(sqrt(post_weights[[i]][,k])*Z.comb[[i]]$z.comb)
      })

      #browser()
      S.tilde.k = Reduce("+",S.tilde.k)/(N*pi[k])

      # S.tilde_k = 0
      # for(i in 1:N){
      #   Z.comb.i = Z.comb[[i]]
      #   post_weights.i = post_weights[[i]]
      #   Ji = nrow(Z.comb.i$z.comb)
      #   for(j in 1:Ji){
      #     S.tilde_k = S.tilde_k + post_weights.i[j,k]*tcrossprod(Z.comb.i$z.comb[j,])
      #   }
      # }
      # temp = eigen(S.tilde_k/(N*pi[k]))

      temp = eigen(S.tilde.k)
      Ulist[[k]] = temp$vectors%*%diag(pmax(temp$values-1,0))%*%t(temp$vectors)
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
  pi.est = pi[(usepointmass+length(U.canon)+1):K]

  if(!is.null(adjust)){
    if(adjust == 'prior'){
      if(is.null(nu)){
        nu = R+1
      }
      U.est.adj = lapply(1:length(U.est),function(k){
        nk = N*pi.est[k]
        if(nk>R){
          s2.hat = R*nu/(nk+nu)/sum(diag(solve(nk*(U.est[[k]]+I_R))))
          (nk/(nk+nu-R-1))*U.est[[k]] + ((R+1-nu+s2.hat)/(nk+nu-R-1))*I_R
        }else{
          NULL
        }
      })

    }

    if(adjust == "lb"){
      U.est.adj = lapply(1:length(U.est),function(k){
        nk = N*pi.est[k]
        if(nk>1){
          U.est[[k]] + 2*diag(sqrt(1/(2*nk)),R)
          #U.est[[k]] + 2*diag(sqrt(2/nk*1/(diag(solve(U.est[[k]]+diag(R))))^2))
          #temp = eigen(U.est[[k]])
          #temp$vectors%*%(diag(pmax(2/sqrt(nk),temp$values)))%*%t(temp$vectors)
        }else{
          NULL
        }
      })
    }
  }

  U.est.adj = U.est.adj[lengths(U.est.adj) != 0]

  # add eb step to adjust estimated U


  result = list()
  result$U.est = U.est
  result$U.est.adj = U.est.adj
  result$pi = pi
  result$loglik = loglik.obj

  result
}
