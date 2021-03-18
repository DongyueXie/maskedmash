

#'@title masked mash for FDR control
#'@description order the tests from most to least significant based on lfsr, and reject least significant ones until FDP<alpha
#'@param Z matrix of Z scores
#'@param obj fitted masked mash object
#'@param alpha target FDR level

maskedmashFDR = function(Z,obj,alpha = 0.05,P=NULL){


  if(is.null(P)){
    P = 2*(1-pnorm(abs(Z)))
  }

  N = nrow(P)
  R = ncol(P)


  lfsr = obj$result$lfsr

  # rank the tests from most significant to least significant
  # and also as initial rejection set
  rej.set = order(lfsr)


  fdp.t = fdp.hat(Z,rej.set,P)
  fdp.tr = fdp.t
  rej.idx = length(rej.set)
  rej.idx.lb = 1
  rej.idx.rb = length(rej.set)

  if(fdp.t>alpha){

    while(fdp.t>alpha | (fdp.t<=alpha&fdp.tr<=alpha)){

      if(rej.idx.lb==rej.idx.rb){
        rej.idx = NULL
        break
      }

      rej.idx = floor((rej.idx.lb+rej.idx.rb)/2)
      fdp.t = fdp.hat(Z,rej.set[1:rej.idx],P)
      fdp.tr = fdp.hat(Z,rej.set[1:(rej.idx+1)],P)

      if(fdp.t>alpha){
        rej.idx.rb = rej.idx
      }

      if((fdp.t<=alpha&fdp.tr<=alpha)){
        rej.idx.lb = rej.idx
      }

    }

  }

  if(is.null(rej.idx)){
    list(rej.set=NULL,fdp=0)
  }else{
    list(rej.set=rej.set[1:rej.idx],fdp=fdp.t)
  }

}


fdp.hat = function(Z,rej.set,P){

  if(is.null(P)){
    P = 2*(1-pnorm(abs(Z)))
  }

  R = sum(P[rej.set]<=1/2)
  A = sum(P[rej.set]>1/2)


  (1+A)/(max(R,1))

}
