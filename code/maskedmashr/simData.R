
#'@titlde simulate from a mash model
#'@description If V=NULL, use V = error_sd^2 I; otherwise, use input as V.
#'@param N total number of samples
#'@param Ulist list of correlation matrices
#'@param pi length K vector
#'@param error_sd standard error of error term.
#'@param signal_sd standard error of error term.
#'@param V array of correlation matrices
#'@return a list of: B, true effect; Bhat: true effect + noise.

simData = function(N,Ulist,V,pi=NULL,error_sd = 1,signal_sd = 3){

  K = length(Ulist)
  R = ncol(Ulist[[1]])
  if(is.null(pi)){pi = rep(1/K,K)}
  nsamp = round(N*pi)
  B = c()
  Bhat = c()
  for(k in 1:K){
    Bk = rmvnorm(nsamp[k],rep(0,R),signal_sd^2*Ulist[[k]])
    B = rbind(B,Bk)
  }

  N = sum(nsamp)
  Bhat = matrix(nrow=N,ncol=R)
  for(i in 1:N){
    Bhat[i,] = rmvnorm(1,B[i,],error_sd^2*V[,,i])
  }
  Shat = matrix(error_sd,nrow=sum(nsamp),ncol=R)

  list(B=B,Bhat=Bhat,Shat=Shat,V=V)

}
