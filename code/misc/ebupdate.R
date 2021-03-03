ebupdate = function(Z.comb,U.canon,U.data,grid,usepointmass,pi,max_iter,tol,
                    U.update,pi_thresh,normalizeU,verbose,printevery){

  N = length(Z.comb)
  R = ncol(U.data[[1]])
  grid.ratio = grid/(1+grid)
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

  ### all possible z-score combinations

  # Z.list = lapply(seq_len(nrow(Z)),function(i){Z[i,]})
  # Z.comb = lapply(Z.list, function(x){enumerate.z(x,thresh)})


  # calculate marginal likelihood
  lik.list = lapply(Z.comb,function(x){
    temp = calc_relative_lik_matrix(mash_set_data(x$z.comb),xUlist,algorithm.version = 'Rcpp')
    temp
  })

  loglik.obj = -Inf

  for(iter in 1:max_iter){

    if(U.update=='rank1'){
      # calculate marginal likelihood
      lik.list = lapply(Z.comb,function(x){
        temp = calc_relative_lik_matrix(mash_set_data(x$z.comb),xUlist,algorithm.version = 'Rcpp')
        temp
      })
    }




    # evaluate objective function
    loglik.obj[iter+1] = sum(log(unlist(lapply(1:N,function(i){
      lm_i = lik.list[[i]]
      det_i = Z.comb[[i]]$z.det.jacob
      sum(det_i*exp(compute_vloglik_from_matrix_and_pi(pi,lm_i, matrix(1,nrow=length(det_i),ncol=R))))
    }))))


    ## check convergence
    if((loglik.obj[iter+1]-loglik.obj[iter])<tol){
      break
    }

    # update gamma_ijkl
    post_weights = lapply(1:N,function(i){
      lik.mat = exp(lik.list[[i]]$loglik_matrix+lik.list[[i]]$lfactors)
      gamma_i = t(t(lik.mat*Z.comb[[i]]$z.det.jacob)*pi)
      gamma_i/sum(gamma_i)
    })


    ## update pi
    pi = lapply(post_weights,colSums)
    pi = Reduce('+',pi)/N

    # update U.data
    ## eigenvectors of sum_ij(S_ij sum_l gamma_ijkl w_l/(1+w_l))

    ## let gamma.tilde_ij = sum_l gamma_ijkl w_l/(1+w_l)

    if(U.update=='rank1'){

      for(k in (length(U.canon)+1):K){

        # find k_index: index the position of k in post_weigfhts matrix
        k_index = usepointmass+seq(k,(n_grid*K),by = K)
        gamma.tilde_k = lapply(post_weights,function(x){
          # figure out gamma_kl
          colSums(t(x[,k_index])*grid.ratio)
        })

        S.tilde_k = 0
        for(i in 1:N){
          Ji = nrow(Z.comb[[i]]$z.comb)
          for(j in 1:Ji){
            S.tilde_k = S.tilde_k + gamma.tilde_k[[i]][j]*tcrossprod(Z.comb[[i]]$z.comb[j,])
          }
        }

        uk = eigen(S.tilde_k)$vectors[,1]
        Ulist[[k]] = tcrossprod(uk)

      }

      ## update xUlist
      xUlist = expand_cov(Ulist,grid,usepointmass)

    }

    if(verbose){
      if(iter%%printevery==0){
        print(sprintf("Done %1.0f iterations, loglikelihood at %.2f",iter,loglik.obj[iter+1]))
      }
    }

  }

  result = calc_post_summary(Z.comb,xUlist,pi,post_weights,pi_thresh)
  result$PosteriorWeights = post_weights
  result$PriorCov = Ulist
  result$pi = pi
  result$grid = grid
  result$loglik = loglik.obj

  result

}
