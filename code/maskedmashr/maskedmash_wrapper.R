

#'@param Bhat matrix of test statistics
#'@param P matrix of p values

maskedmash_wrapper = function(Bhat,P=NULL,thresh=NULL,npc=5,Shat=NULL,adjust = 'lb',
                              nu = ncol(Bhat)+1,verbose=FALSE){

  R = ncol(Bhat)
  N = nrow(Bhat)
  if(!is.null(P)){
    Z = sign(Bhat)*qnorm(1-P/2)
  }else{
    Z = Bhat
  }
  data = mash_set_data(Z,Shat)
  U.c = cov_canonical(data)
  m.1by1 = mash_1by1(data)
  strong = get_significant_results(m.1by1)
  if(length(strong)==0){
    strong=NULL
  }
  U.pca = cov_pca(data,npc,strong)
  #browser()
  U.est = masked.md(data,strong=strong,thresh=thresh,U.data=U.pca,adjust=adjust,verbose=verbose)$U.est.adj
  out = masked.mash(data,thresh=thresh,U.canon = U.c,U.data = U.est,verbose=verbose)
  out

}
