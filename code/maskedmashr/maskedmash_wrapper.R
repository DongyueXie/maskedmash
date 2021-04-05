


#'@description
#'@param Bhat matrix of estimated effects/Z scores
#'@param P matrix of p values
#'@param p.thresh threshold of p-values, take values in [0,0.5]. p<=thresh and p>=(1-thresh) will be masked.

maskedmash_wrapper = function(Bhat,P,p.thresh=0.5,npc=5,Shat=NULL,adjust = 'lb',verbose=FALSE,eps=1e-15){

  R = ncol(Bhat)
  N = nrow(Bhat)

  # if(!is.null(P)){
  #   Z = sign(Bhat)*qnorm(1-P/2)
  # }else{
  #   Z = Bhat
  # }

  P.temp = P
  P.temp = pmax(P.temp,eps)
  P.temp = pmin(P.temp,1-eps)
  #P.temp[which(P.temp<eps)] = eps
  #P.temp[which(P.temp>(1-eps))] = 1-eps

  Z = sign(Bhat)*qnorm(1-P.temp/2)
  if(is.null(p.thresh)){
    p.thresh = 0.5
  }
  thresh = qnorm(1-p.thresh/2)


  data = mash_set_data(Z,Shat)
  U.c = cov_canonical(data)
  m.1by1 = mash_1by1(data)
  strong = get_significant_results(m.1by1,0.2)
  if(length(strong)==0){
    strong=NULL
  }
  U.pca = cov_pca(data,npc,strong)
  #browser()
  U.est = masked.md(data,strong=strong,thresh=thresh,U.data=U.pca,adjust=adjust,verbose=verbose)$U.est.adj
  out = masked.mash(data,thresh=thresh,U.canon = U.c,U.data = U.est,verbose=verbose)
  out$p.thresh = p.thresh
  out$P = P
  out
}
