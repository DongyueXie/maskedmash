
estimate_pi = function(Z.comb,U.canon,U.data,grid,usepointmass,pi,max_iter,tol,
                    pi_thresh,normalizeU,algorithm.version,optmethod,verbose,printevery,control,prior){

  N = length(Z.comb)
  R = ncol(U.data[[1]])
  n_grid = length(grid)

  Ulist = c(U.canon,U.data)
  if(normalizeU){Ulist = normalize_Ulist(Ulist)}
  K = length(Ulist)
  xUlist = expand_cov(Ulist,grid,usepointmass)
  P = length(xUlist) # the order of P is : null, (l=1,k=1:K), (l=2,k=1:K),...
  if(is.null(pi)){pi = rep(1/P,P)}
  # marginal likelihood
  ## need to include L
  ### marginal likelihood list: a list of length N, each element is a matrix of dim Ji*P.
  ### posterior weights: a list of length N, each element is a matrix of dim Ji*P
  ### prior weights pi is a length P vector

  # calculate likelihood
  lik.list = lapply(Z.comb,function(x){
    temp = calc_relative_lik_matrix(mash_set_data(x$z.comb),xUlist,algorithm.version = algorithm.version)
    temp
  })

  if(optmethod=='mixSQP'){

    # Get the likelihood matrix, of dimension (sum_i Ji) * K
    lik_mat = lapply(1:N,function(i){
      colSums(Z.comb[[i]]$z.det.jacob * exp(lik.list[[i]]$loglik_matrix+lik.list[[i]]$lfactors))
      })
    lik_mat = do.call(rbind,lik_mat)
    prior = set_prior(ncol(lik_mat),prior)
    pi = optimize_pi(lik_mat,prior=prior,optmethod=optmethod,control=control)

  }

  if(optmethod=='EM'){

    if(prior == 'uniform'){
      multiplier = rep(1,P)
    }
    if(prior == 'nullbiased'){
      multiplier = c(10,rep(1,P-1))
    }

    loglik.obj = -Inf

    for(iter in 1:max_iter){

      # evaluate objective function
      loglik.obj[iter+1] = calc_obj(lik.list,Z.comb,pi,N,R,multiplier)

      ## check convergence
      if((loglik.obj[iter+1]-loglik.obj[iter])<tol){
        break
      }

      # update gamma_ijkl
      post_weights = calc_post_weights(Z.comb,lik.list,pi)


      ## update pi
      pi = lapply(post_weights,colSums)
      nk = Reduce('+',pi) + multiplier - 1
      pi = nk/sum(nk)

      if(verbose){
        if(iter%%printevery==0){
          print(sprintf("Done %1.0f iterations, loglikelihood at %.2f",iter,loglik.obj[iter+1]))
        }
      }

    }
  }



  #result = calc_post_summary(Z.comb,xUlist,pi,post_weights,pi_thresh)
  result = list()
  result$post_weights = calc_post_weights(Z.comb,lik.list,pi)
  result$xUlist = xUlist
  result$pi = pi
  result$loglik = calc_obj(lik.list,Z.comb,pi,N,R)

  result

}
