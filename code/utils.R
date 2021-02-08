
#'@title masking function
mask.func = function(z){
  qnorm(1.5-pnorm((z)))
}

#'@title returns all possible combinations of masked z-scores and |det(Jacobian)|
#'@param z a vector
#'@param thresh |z| larger than thresh or s(|z|)<=s(t) will be masked. when thresh<=Phi^{-1}(0.75), all masked.
#'@return a J*R matrix of all possible z scores and a length J vector of |det(Jacobian)|.
enumerate.z = function(z,thresh){
  is.mask = (abs(z)>=thresh)|(abs(z)<=mask.func(thresh))
  s_z = sign(z)*mask.func(abs(z))*is.mask + z*(1-is.mask)
  z.mask = pmin(z,s_z)
  all.comb = as.matrix(expand.grid(as.list(as.data.frame(rbind(z,s_z)))))
  z.comb = all.comb[!duplicated(all.comb),,drop=FALSE]
  z.det.jacob = apply(z.comb,1,function(x){abs(prod(dnorm(z.mask)/dnorm(x)))})
  return(list(z.mask=z.mask,z.comb = z.comb,z.det.jacob=z.det.jacob))
}



#'@title set augmented z scores
#'@param Z z score matrix
#'@param thresh
#'@reutrn a list, Z.comb: each element is a Ji*R matrix; z.mask: masked z scores; z.det.jacob
set_Z.comb = function(Z,thresh){
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
  list(fdp=ff,power=pp)
}

rmse = function(x,y){sqrt(mean((x-y)^2))}
