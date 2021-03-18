
#'@param N total number of samples
#'@param Ulist a list of correlation matrices
#'@param pi prior weights
#'@param prior normal, t, uniform
#'@param signal_sd variance of mu's, used as signal_sd^2*U
#'@param half.uniform True to generate only positive uniforms
#'@param unif.range range of uniform priors
simDataI.ult = function(N,Ulist,pi=NULL,prior="t",
                        signal_sd = sqrt(3),
                        df=10,half.uniform=FALSE,
                        unif.range = 3){


  K = length(Ulist)
  R = ncol(Ulist[[1]])
  if(is.null(pi)){pi = rep(1/K,K)}
  nsamp = round(N*pi)
  B = c()
  Bhat = c()

  if(prior=='normal'){

    for(k in 1:K){
      Bk = rmvnorm(nsamp[k],rep(0,R),signal_sd^2*Ulist[[k]])
      B = rbind(B,Bk)
    }

  }else if(prior == 't'){

    for(k in 1:K){
      if(df>2){
        Bk = mvtnorm::rmvt(nsamp[k],sigma=signal_sd^2*Ulist[[k]],df=df)
      }else{
        Bk = mvtnorm::rmvt(nsamp[k],sigma=signal_sd^2*Ulist[[k]],df=df)
      }

      B = rbind(B,Bk)
    }

  }else if(prior == 'uniform'){

    for(k in 1:K){
      Bk = rmvu(nsamp[k],Ulist[[k]])
      B = rbind(B,Bk)
    }

    #temp = sqrt(signal_sd^2*12)

    if(half.uniform){
      B = B*unif.range
    }else{
      non0.idx = which(B!=0)
      B[non0.idx] = B[non0.idx]*2*unif.range-unif.range
    }


  }else{
    stop('prior should be one of normal, t, uniform')
  }


  N = sum(nsamp)
  # Bhat = matrix(nrow=N,ncol=R)
  # for(i in 1:N){
  #   Bhat[i,] = rmvnorm(1,B[i,],diag(R))
  # }
  Bhat = matrix(rnorm(N*R,B,sd=1),ncol = R)
  Shat = matrix(1,nrow=N,ncol=R)

  list(B=B,Bhat=Bhat,Shat=Shat)


}




rmvu = function(n,U){
  non0.idx = which(diag(U)!=0)
  if(length(non0.idx)==ncol(U)){
    U = cov2cor(U)
    X = mvtnorm::rmvnorm(n,rep(0,ncol(U)),U)
    X = apply(X,2,pnorm)
  }else if(length(non0.idx)==0){
    X = matrix(0,nrow=n,ncol=ncol(U))
  }else{
    U.non0 = U[non0.idx,non0.idx,drop=F]
    U.non0 = cov2cor(U.non0)
    X = matrix(0,nrow=n,ncol=ncol(U))
    x = mvtnorm::rmvnorm(n,rep(0,ncol(U.non0)),U.non0)
    x = apply(x,2,pnorm)
    X[,non0.idx] = x
  }
  X
}



