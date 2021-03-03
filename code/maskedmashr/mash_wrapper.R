

mash_wrapper = function(Bhat,P=NULL,npc=5,Shat=NULL,nu = ncol(Bhat)+1,verbose=FALSE){



  if(!is.null(P)){
    Z = sign(Bhat)*qnorm(1-P/2)
  }else{
    Z = Bhat
  }

  R = ncol(Z)
  N = nrow(Z)

  data = mash_set_data(Z,Shat)
  U.c = cov_canonical(data)
  m.1by1 = mash_1by1(data)
  strong = get_significant_results(m.1by1)
  n.ed = length(strong)
  if(length(strong)==0){
    strong=NULL
    n.ed=N
  }
  U.pca = cov_pca(data,npc,strong)
  U.ed = bovy_wrapper(data,U.pca,subset=strong)
  U.est = U.ed$Ulist
  U.est = lapply(1:length(U.est),function(k){
    nk = n.ed*U.ed$pi[k]
    if(nk>R){
      s2.hat = R*nu/(nk+nu)/sum(diag(solve(nk*(U.est[[k]]+diag(R)))))
      (nk/(nk+nu-R-1))*U.est[[k]] + ((R+1-nu+s2.hat)/(nk+nu-R-1))*diag(R)
    }else{
      NULL
    }

  })

  U.est = U.est[lengths(U.est) != 0]


  out = mash(data,c(U.c,U.est),verbose=verbose)
  out

}
