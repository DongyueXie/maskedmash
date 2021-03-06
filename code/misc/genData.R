
#'@param n number of samples for each cov
#'@param Ulist covriance matrices

genData = function(n,Ulist,signal_sd,error_sd=1){
  R = nrow(Ulist[[1]])
  K = length(Ulist)
  B = c()
  Bhat = c()
  for(k in 1:K){
    Bk = rmvnorm(n,sigma = signal_sd^2*Ulist[[k]])
    Xk = t(apply(Bk,1,function(x){rmvnorm(1,x,error_sd^2*diag(R))}))
    B = rbind(B,Bk)
    Bhat = rbind(Bhat,Xk)
  }

  return(list(B = B, Bhat= Bhat,Shat = matrix(error_sd,nrow=K*n,ncol=R)))
}
