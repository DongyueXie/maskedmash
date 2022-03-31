args = commandArgs(trailingOnly=TRUE)
# determine the value of epsilon to add to the diagonal of prior covariances
epsilon = as.numeric(args[1])

### implement mash with gene-wise true median subtracted from bhat
setwd("/scratch/midway2/yushaliu/multivariate_poisson/simulations_sc/setting2/signal_mash")
library(Matrix)
library(mashr)

# load in the true prior covariance
ulist <- readRDS("sim_design.rds")$ulist


for(idx in 1:10){
  ################################ Read in data #########################################
  data <- readRDS(file=paste0("data/iter", idx, ".Rds"))
  bhat <- data$bhat
  shat <- data$shat
  R <- ncol(bhat)
  bhat <- bhat - data$mu %*% t(rep(1,R))

  
  ################################ Run mash by directly subtracting the median ##################################
  start_time = proc.time()
  print("##########################################")
  print(sprintf("start mash with true median subtracted directly for replicate %d", idx))
  
  ### set up for mash 
  mash.data = mash_set_data(bhat, shat, alpha = 1)
  
  ### use true prior covariances
  U = list()
  U[[1]] <- ulist[[1]] %*% t(ulist[[1]])
  U[[2]] <- ulist[[2]] %*% t(ulist[[2]])
  U[[3]] <- ulist[[3]] %*% t(ulist[[3]])
  names(U) <- c("PC1", "PC2", "PC3")
  for(k in 1:length(U)){
    Uk <- U[[k]]
    U[[k]] <- max(diag(Uk))*(Uk/max(diag(Uk)) + epsilon*diag(ncol(Uk)))
  }
  
  ### fit mash
  mash.fit = mash(mash.data, U, algorithm.version = 'R')
  print(sprintf("finish mash with true median subtracted directly for replicate %d", idx))
  runtime = proc.time() - start_time
  mash.fit[["runtime"]] = runtime
  saveRDS(mash.fit, file = paste0("output/mash_oracle_epsilon_", epsilon, "_rep_", idx, ".Rds"))  
  
  rm(data, bhat, shat, mash.fit)
}



print(sessionInfo())