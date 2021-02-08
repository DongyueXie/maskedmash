
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
#'@param verbose TRUE to print progress
#'@return a list of Posterior mean, sd, lfsr, lfdr, negativeProb.
masked.mash = function(data,thresh=NULL,pi=NULL,U.canon=NULL,U.data,
                       fixg=TRUE,usepointmass=TRUE,
                       U.update = 'none',
                       pi_thresh = 1e-5,
                       gridmult = sqrt(2),
                       grid = NULL,
                       normalizeU = TRUE,
                       max_iter=1000,
                       tol=1e-3,
                       verbose=TRUE,
                       printevery = 1){
  Z = data$Bhat
  N = nrow(Z)
  R = ncol(Z)
  I_R = diag(R)

  if(is.null(thresh)){
    # mask all a-scores
    thresh = qnorm(0.75)
  }
  Z.comb = set_Z.comb(Z,thresh)

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
      temp = calc_relative_lik_matrix(mash_set_data(x$z.comb),xUlist)
      temp
    })

    # calculate posterior weights
    post_weights = calc_post_weights(Z.comb,lik.list,pi)
    #browser()
    # calculate posteriors
    result = calc_post_summary(Z.comb,xUlist,pi,post_weights,pi_thresh)
    result$PriorCov = Ulist
    result$PosteriorWeights = post_weights

  }else{

    ############## estimate weights and covariances ####################
    if(is.null(grid)){
      grid = autoselect_grid(data,gridmult)
    }
    result = ebupdate(Z.comb,U.canon,U.data,grid,usepointmass,
                      pi,max_iter,tol,U.update,pi_thresh,normalizeU,verbose,printevery)


  }


  effect_names = rownames(data$Bhat)
  condition_names = colnames(data$Bhat)
  for (i in names(result)) {
    if (length(dim(result[[i]])) == 2) {
      colnames(result[[i]]) = condition_names
      rownames(result[[i]]) = effect_names
    }
  }

  result$maskedProp = sum(log(unlist(lapply(Z.comb,function(x){nrow(x$z.comb)})),2))/(N*R)


  result


}







#
# masked.mash = function(X,w,U){
#   N = nrow(X)
#   R = ncol(X)
#   K = length(U)
#   S = diag(R)
#   result = vector(mode='list',length = 6)
#   names(result) = c('PosteriorMean','PosteriorSD','lfdr','lfsr','NegativeProb','PosteriorWeight')
#   result = lapply(result,function(z){matrix(nrow=N,ncol=R)})
#   result$PosteriorWeight = matrix(nrow=N,ncol=2^R*K)
#
#   xUlist = expand_cov(U,1,FALSE)
#
#   for(i in 1:N){
#     # for each sample, generate all possible 'observations'
#     xi = enumerate.recip(X[i,])
#
#     ## 1. calculate each observation's likelihood hence posterior weights
#     datai = mash_set_data(xi)
#     lmi = calc_relative_lik_matrix(datai,xUlist)
#     xi_llik = compute_vloglik_from_matrix_and_pi(w,lmi,datai$Shat_alpha)
#     w_xi = exp(xi_llik-max(xi_llik))
#     w_xi = c(w_xi/sum(w_xi))
#     ## 2. then conditional on each likelihood, obtain posteriors
#     w_mixi = compute_posterior_weights(w,exp(lmi$loglik_matrix))
#     mat_posti = compute_posterior_matrices(datai, xUlist,w_mixi)
#     ### now we have weights for observations and each mixture component.
#     mat_posti = lapply(mat_posti,function(z){colSums(z*w_xi)})
#     result$PosteriorMean[i,] = mat_posti$PosteriorMean
#     result$PosteriorSD[i,] = mat_posti$PosteriorSD
#     result$lfdr[i,] = mat_posti$lfdr
#     result$lfsr[i,] = mat_posti$lfsr
#     result$NegativeProb[i,] = mat_posti$NegativeProb
#   }
#
#   result
#
# }
#
# enumerate.sign = function(x){
#   x.sign = rbind(x,-x)
#   as.matrix(expand.grid(as.list(as.data.frame(x.sign))))
# }
#
# enumerate.recip = function(x){
#   x.sign = rbind(x,1/x)
#   as.matrix(expand.grid(as.list(as.data.frame(x.sign))))
# }
#
#
#
# llik_tt = c()
# for(j in 1:32){
#   llik_tt[j] = (dmvnorm(xi[j,],sigma=diag(5)+U[[1]]) + dmvnorm(xi[j,],sigma=diag(5)+U[[2]]))/2
# }
# llik_tt/sum(llik_tt)
