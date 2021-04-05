
#'@title masking function
mask.func = function(z){
  sign(z)*qnorm(1.5-pnorm(abs(z)))
}

#'@title returns all possible combinations of masked z-scores and |det(Jacobian)|
#'@param z a vector
#'@param thresh |z| larger than thresh or s(|z|)<=s(t) will be masked. when thresh<=Phi^{-1}(0.75), all masked.
#'@return a J*R matrix of all possible z scores, a length J vector of |det(Jacobian)|, and masked z scores.
enumerate.z = function(z,thresh){
  masking = is.mask.z(z,thresh)
  s_z = mask.func(z)*masking + z*(1-masking)
  z.mask = pmin(z,s_z)
  all.comb = as.matrix(expand.grid(as.list(as.data.frame(rbind(z,s_z)))))
  z.comb = all.comb[!duplicated(all.comb),,drop=FALSE]
  z.det.jacob = apply(z.comb,1,function(x){
    exp(sum(dnorm(z.mask,log=TRUE) - dnorm(x,log=TRUE)))
    #abs(prod(dnorm(z.mask[dif.idx])/dnorm(x[dif.idx])))
    })
  return(list(z.mask=z.mask,z.comb = z.comb,z.det.jacob=z.det.jacob))
}

is.mask.z = function(z,thresh){
  ((abs(z)>=thresh)|(abs(z)<=mask.func(thresh)))
}

# is.mask.z = function(z,thresh,max_mask = 10){
#   ((abs(z)>=thresh)|(abs(z)<=mask.func(thresh)))&((abs(z)<max_mask)&(abs(z)>mask.func(max_mask)))
# }

is.mask.p = function(p,thresh){
  if(thresh<=0.5){
    (p<=thresh)|(p>=(1-thresh))
  }else{
    stop('thresh should in [0,0.5]')
  }

}



#'@title set augmented z scores
#'@param Z z score matrix
#'@param thresh
#'@reutrn a list, Z.comb: each element is a Ji*R matrix; z.mask: masked z scores; z.det.jacob
mask.Z = function(Z,thresh){
  Z.list = lapply(seq_len(nrow(Z)),function(i){Z[i,]})
  Z.comb = lapply(Z.list, function(x){enumerate.z(x,thresh)})
  Z.comb
}



fdp = function(dis.idx, true.idx){
  if(length(dis.idx)==0){
    0
  }else{
    1-mean(dis.idx%in%true.idx)
  }
}


#'@param bhat rejected statisitcs
#'@param b true b of those rejected ones.
fsp = function(bhat,b){
  if(length(bhat)==0){
    0
  }else{
    sum(sign(bhat)!=sign(b))/length(bhat)
  }
}


auc = function(pred,true.label){
  auc=pROC::roc(response = true.label, predictor = pred,direction = '<',levels = c(0,1))
  auc$auc
}

powr = function(dis.idx, true.idx){
  if(length(dis.idx)==0){
    0
  }else{
    sum(dis.idx%in%true.idx)/length(true.idx)
  }
}

fdp_power = function(lfsr,B,alpha=0.05){
  non_null_idx = which(B!=0)
  ff = fdp(which(lfsr<alpha),non_null_idx)
  pp = powr(which(lfsr<alpha),non_null_idx)
  #fdp.hat = mean(lfsr[which(lfsr<alpha)])
  list(fdp=ff,power=pp)
}


fsp_power = function(lfsr,Bhat,B,alpha=0.05){
  non_null_idx = which(B!=0)
  dis.idx = which(lfsr<alpha)
  ff = fsp(Bhat[dis.idx],B[dis.idx])
  pp = powr(dis.idx,non_null_idx)
  fsp.hat = mean(lfsr[dis.idx])
  list(fsp=ff,power=pp,fsp.hat=fsp.hat)
}

rmse = function(x,y){sqrt(mean((x-y)^2))}
