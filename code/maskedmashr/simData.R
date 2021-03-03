
#'@titlde simulate from a mash model
#'@param N total number of samples
#'@param Ulist list of covariance matrices
#'@param pi length K vector
#'@param error_sd standard error of error term.
#'@return a list of: B, true effect; Bhat: true effect + noise.

simData = function(N,Ulist,pi=NULL,err_sd = 1){

  K = length(Ulist)
  R = ncol(Ulist[[1]])
  if(is.null(pi)){pi = rep(1/K,K)}
  nsamp = round(N*pi)
  B = c()
  Bhat = c()
  for(k in 1:K){
    Bk = rmvnorm(nsamp[k],rep(0,R),Ulist[[k]])
    B = rbind(B,Bk)
  }

  N = sum(nsamp)
  Bhat = matrix(nrow=N,ncol=R)
  for(i in 1:N){
    Bhat[i,] = rmvnorm(1,B[i,],err_sd^2*diag(R))
  }

  list(B=B,Bhat=Bhat,Shat = matrix(err_sd,nrow=sum(nsamp),ncol=R))

}
