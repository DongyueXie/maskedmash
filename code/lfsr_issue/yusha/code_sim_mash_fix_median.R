args = commandArgs(trailingOnly=TRUE)
iter.idx = as.integer(args[1])

### implement mash with gene-wise median estimated and subtracted from bhat
setwd("/scratch/midway2/yushaliu/multivariate_poisson/simulations_sc/setting2/signal_mash")
library(Matrix)
library(mashr)

# determine the value of epsilon to add to the diagonal of prior covariances
epsilon <- 1e-2


for(idx in iter.idx:(iter.idx+1)){
  ################################ Read in data #########################################
  data <- readRDS(file=paste0("data/iter", idx, ".Rds"))
  bhat <- data$bhat
  shat <- data$shat
  R <- ncol(bhat)
  bhat <- bhat - apply(bhat, 1, median) %*% t(rep(1,R))

  
  ################################ Run mash by directly subtracting the median ##################################
  start_time = proc.time()
  print("##########################################")
  print(sprintf("start mash with median estimated and subtracted directly for replicate %d", idx))
  
  ### set up for mash 
  mash.data = mash_set_data(bhat, shat, alpha = 1)
  
  ### canonical covariances
  U.c = cov_canonical(mash.data)
  for(k in 1:length(U.c)){
    Uk <- U.c[[k]]
    U.c[[k]] <- max(diag(Uk))*(Uk/max(diag(Uk)) + epsilon*diag(ncol(Uk)))
  }
  
  ### data-driven covariances
  m.1by1 = mash_1by1(mash.data, alpha=1)
  strong = get_significant_results(m.1by1, 0.05)
  U.pca = cov_pca(mash.data, 5, subset=strong)
  U.ed = cov_ed(mash.data, U.pca, subset=strong) 
  for(k in 1:length(U.ed)){
    Uk <- U.ed[[k]]
    U.ed[[k]] <- max(diag(Uk))*(Uk/max(diag(Uk)) + epsilon*diag(ncol(Uk)))
  }
  
  ### fit mash
  mash.fit = mash(mash.data, c(U.c,U.ed), algorithm.version = 'R')
  print(sprintf("finish mash with median estimated and subtracted directly for replicate %d", idx))
  runtime = proc.time() - start_time
  mash.fit[["runtime"]] = runtime
  saveRDS(mash.fit, file = paste0("output/mash_fix_median_rep_", idx, ".Rds"))  
  saveRDS(list(U.c=U.c, U.ed=U.ed), file = paste0("output/mash_fix_median_U_rep_", idx, ".Rds"))  
  
  rm(data, bhat, shat, mash.fit)
}



print(sessionInfo())