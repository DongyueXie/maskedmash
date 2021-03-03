

#'@title mask on masked data
#'@description The procedure masks the z scores, by taking the bigger one(absolute value) in masked set, then fit mash.


mash_mask = function(Bhat,P=NULL,thresh=NULL,npc=5,Shat=NULL,nu = ncol(Z)+1,verbose=FALSE){

  if(is.null(thresh)){
    thresh = 0
  }
  if(!is.null(P)){
    Z = sign(Bhat)*qnorm(1-P/2)
  }else{
    Z = Bhat
  }
  N = nrow(Z)
  R = ncol(Z)
  # mask z scores
  Z.mask = t(apply(Z,1,function(x){
    x.alt = sign(x)*qnorm(1.5-pnorm(abs(x)))
    is.mask = (abs(x)>=thresh)|(abs(x)<=mask.func(thresh))
    is.mask*sign(x)*pmax(abs(x),abs(x.alt)) + (1-is.mask)*x
  }))

  mash_wrapper(Bhat=Z.mask,npc=npc,Shat=Shat,nu=nu,verbose=verbose)


}
