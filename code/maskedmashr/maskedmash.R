
## implementation ##

#'@title masked mash
#'@param data object from mash_set_data
#'@param thresh NULL or a number, |z-score| larger than thresh will be masked.
#'@param pi prior weights; either fixed or as init value
#'@param U.canon a list of canonical prior cov matrices; fixed.
#'@param U.data a list of data-driven prior cov marrices; either fixed or as initialization.
#'@param fixg whether fix the prior weights and covariance matrices.
#'@param usepointmass whether to include a point mass at 0
#'@param pi_thresh threshold below which mixture components are ignored
#'@param normalizeU whether normalize U such that the largest diag element is 1
#'@param algorithm.version Rcpp or R for evaluating likelihood
#'@param optmethod 'mixSQP' or 'EM'
#'@param verbose TRUE to print progress
#'@param control a list of control parameters for SQP
#'@param prior nullbiased or uniform
#'@return a list of Posterior mean, sd, lfsr, lfdr, negativeProb.
masked.mash = function(data,thresh=NULL,pi=NULL,U.canon=NULL,U.data,
                       fixg=FALSE,usepointmass=TRUE,
                       #U.update = 'none',
                       pi_thresh = 1e-5,
                       gridmult = sqrt(2),
                       grid = NULL,
                       normalizeU = TRUE,
                       algorithm.version = 'Rcpp',
                       optmethod = 'mixSQP',
                       prior = 'nullbiased',
                       max_iter=1000,
                       tol=1e-3,
                       verbose=TRUE,
                       printevery = 100,
                       control=list()){
  Z = data$Bhat
  N = nrow(Z)
  R = ncol(Z)
  I_R = diag(R)

  if(is.null(thresh)){
    # mask all a-scores
    thresh = qnorm(0.75)
  }
  Z.comb = mask.Z(Z,thresh)

  # if prior is fixed, then directly calculate posterior summaries

  ## posterior weights: of dimension N*J_i*K. Use list, each list element is a J_i*K matrix
  ##
  if(fixg){

    Ulist = c(U.canon,U.data)
    if(normalizeU){Ulist = normalize_Ulist(Ulist)}
    if(!is.null(grid)){
      xUlist = expand_cov(Ulist,grid,usepointmass)
    }else{
      if(usepointmass){
        xUlist = c(list(null=matrix(0,nrow = R,ncol = R)),Ulist)
      }else{
        xUlist = Ulist
      }
    }

    if(is.null(pi)){pi = rep(1/length(xUlist),length(xUlist))}
    lik.list = lapply(Z.comb,function(x){
      temp = calc_relative_lik_matrix(mash_set_data(x$z.comb),xUlist,algorithm.version = algorithm.version)
      temp
    })

    # calculate posterior weights
    post_weights = calc_post_weights(Z.comb,lik.list,pi)
    #browser()
    # calculate posteriors
    #result = calc_post_summary(Z.comb,xUlist,pi,post_weights,pi_thresh)
    #result$PriorCov = Ulist
    #result$PosteriorWeights = post_weights
    loglik = calc_obj(lik.list,Z.comb,pi,N,R)

  }else{

    ############## estimate weights ####################
    if(is.null(grid)){
      grid = autoselect_grid(data,gridmult)
    }
    out = estimate_pi(Z.comb,U.canon,U.data,grid,usepointmass,
                      pi,max_iter,tol,pi_thresh,normalizeU,
                      algorithm.version,optmethod,verbose,printevery,control,prior)
    pi = out$pi
    xUlist = out$xUlist
    post_weights = out$post_weights
    loglik = out$loglik

  }


  # get posterior summaries
  result = calc_post_summary(Z.comb,xUlist,pi,post_weights,pi_thresh)


  effect_names = rownames(data$Bhat)
  condition_names = colnames(data$Bhat)
  for (i in names(result)) {
    if (length(dim(result[[i]])) == 2) {
      colnames(result[[i]]) = condition_names
      rownames(result[[i]]) = effect_names
    }
  }

  out = list()
  out$result = result
  out$maskedProp = sum(log(unlist(lapply(Z.comb,function(x){nrow(x$z.comb)})),2))/(N*R)
  out$loglik = loglik
  out$fitted_g = list(pi=pi,Ulist = c(U.canon,U.data),grid = grid, usepointmass = usepointmass)
  out$posterior_weights = post_weights


  out


}





