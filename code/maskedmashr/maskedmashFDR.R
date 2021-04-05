#'@title masked mash for FDR control
#'@description order the tests from most to least significant based on lfsr, and reject least significant ones until FDP<alpha
#'@param obj fitted masked mash object
#'@param alpha target FDR level

maskedmashFDR = function(obj,alpha = 0.05){

  P = obj$P
  N = nrow(P)
  R = ncol(P)
  score = obj$result$lfsr
  rej.region = which(is.mask.p(P,obj$p.thresh))
  P = P[rej.region]
  score = score[rej.region]

  # rank the tests from most significant to least significant
  # and also as initial rejection set
  rej.set = rej.region[order(score)]


  fdp.t = fdp.hat(rej.set,obj$P)
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
      fdp.t = fdp.hat(rej.set[1:rej.idx],obj$P)
      fdp.tr = fdp.hat(rej.set[1:(rej.idx+1)],obj$P)

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




fdp.hat = function(rej.set,P){

  R = sum(P[rej.set]<=1/2)
  A = sum(P[rej.set]>1/2)
  (1+A)/(max(R,1))

}
























#'@title masked mash for FDR control
#'@description order the tests from most to least significant based on lfsr, and reject least significant ones until FDP<alpha
#'@param Z matrix of Z scores
#'@param obj fitted masked mash object
#'@param alpha target FDR level
#'@param P p-values
#'@param score “score” that measures how ”promising” each hypothesis, smaller score means more likely the hypothesis is non-null.
#'@param rej.region index set of hypothesis in the rejection region(theri p-values being masked).

maskedmashFDR.temp = function(Z = NULL, obj,alpha = 0.05,P=NULL,score = NULL,rej.region = NULL){


  if(is.null(P)){
    if(!is.null(Z)){
      P = 2*(1-pnorm(abs(Z)))
    }else{
      stop("one of P and Z must be given.")
    }
  }

  N = nrow(P)
  R = ncol(P)


  if(is.null(score)){
    if(!is.null(obj)){
      score = obj$result$lfsr
    }else{
      stop("one of score and obj must be given.")
    }
  }

  if(!is.null(rej.region)){
    P = P[rej.region]
    score = score[rej.region]
  }else{
    rej.region = 1:length(P)
  }


  # rank the tests from most significant to least significant
  # and also as initial rejection set
  rej.set = rej.region[order(score)]


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
