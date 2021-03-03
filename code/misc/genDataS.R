genDataS = function(nsamp,U,S,signal_sd){
  K = length(U)
  R = ncol(S)
  B = c()
  Bhat = c()
  for(k in 1:K){
    Bk = mvtnorm::rmvnorm(nsamp,sigma=signal_sd^2*U[[k]])
    B = rbind(B,Bk)
    Bkhat = matrix(nrow=nsamp,ncol=R)
    for(i in 1:nsamp){
      Bkhat[i,] = mvtnorm::rmvnorm(1,mean=Bk[i,],sigma=S)
    }
    Bhat = rbind(Bhat,Bkhat)
  }
  list(B=B,Bhat=Bhat)
}
