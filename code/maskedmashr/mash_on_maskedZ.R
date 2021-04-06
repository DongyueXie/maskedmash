

#'@title mask on masked data
#'@description The procedure masks the z scores, by taking the bigger one(absolute value) in masked set, then fit mash.


mash_mask = function(Bhat,P,p.thresh=0.5,npc=5,Shat=NULL,verbose=FALSE,eps=1e-15){


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

  Z = sign(Bhat)*qnorm(1-P.temp/2)
  if(is.null(p.thresh)){
    p.thresh = 0.5
  }
  thresh = qnorm(1-p.thresh/2)


  # mask z scores
  Z.mask = t(apply(Z,1,function(x){
    x.alt = mask.func(x)
    is.mask = is.mask.z(x,thresh)
    is.mask*sign(x)*pmax(abs(x),abs(x.alt)) + (1-is.mask)*x
  }))

  out = mash_wrapper(Bhat=Z.mask,npc=npc,Shat=Shat,verbose=verbose)
  out$P = P
  out$p.thresh = p.thresh
  out
}
