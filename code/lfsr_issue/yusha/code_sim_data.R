rm(list=ls())
setwd("/Users/mac/Documents/multivariate_poisson/simulations_sc/setting2")
library(Matrix)
library(limma)

################################ Determine simulation parameters based on real data #########################################
scdata <- readRDS("raw_count.rds")
sample.info <- readRDS("sample_info.rds")

### specifiy the selected conditions
trts <- c("CCL11", "CCL17", "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CXCL5", "IFNg", "IL12p70", "IL18",
          "IL1a", "IL1b", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL9", "IL27", "IL34", "IL36a", "TGFb", "TSLP") 

### take the subset of cells from given conditions
idx.cell <- sample.info$sample %in% trts
scdata <- scdata[, idx.cell]
sample.info <- sample.info[idx.cell,]
sum(colnames(scdata)!=sample.info$X0)

### construct the design matrix without intercept
conditions <- factor(droplevels(sample.info$sample), levels=trts)
designmat <- model.matrix(~0+conditions)

### remove genes with totol counts over conditions fewer than 25 before thinning
keep.idx <- which(rowSums(scdata)>=25)
scdata <- scdata[keep.idx,]

### randomly permute the cells with respect to the condition labels independently for each gene
set.seed(1)
scdata <- t(apply(scdata, 1, function(x){x[sample.int(ncol(scdata))]}))

### log transform and normalize single cell count data
scale.factor <- colSums(scdata)
data.mash <- log(t(t(scdata)/scale.factor)*median(scale.factor) + 0.1)

### extract Bhat and Shat using limma  
data.X = model.matrix(~0+conditions)
colnames(data.X) <- trts
cov_of_interest <- 1:ncol(data.X)
lmout <- limma::lmFit(object = data.mash, design = data.X)
eout <- limma::eBayes(lmout)
bhat <- lmout$coefficients[,cov_of_interest,drop=FALSE]
shat <- lmout$stdev.unscaled[,cov_of_interest,drop=FALSE] * sqrt(eout$s2.post)
rownames(shat) <- rownames(bhat)

### determine the gene-wise median
mu <- apply(bhat, 1, median)

### load in the prior covariances which are similar to those of the original data, and do some modifications
ulist.sim <- readRDS("ulist_sim.rds")
ulist <- ulist.sim$ulist
u3 <- ulist[[3]]
u3[u3 < 0.5] <- 0
ulist[[3]] <- u3
pi.u <- rep(1/3, 3)

### save the simulation design parameters
saveRDS(list(mu=mu, shat=shat, ulist=ulist, pi.u=pi.u), "signal_mash/sim_design.rds")



############################## Simulate gaussian data based on the design parameters determined above ###################################
J <- length(mu)
R <- length(trts)

### set the seed
set.seed(100)

for (iter in 1:20){
  ### simulate treatment effects
  beta <- matrix(0, nrow=J, ncol=R)
  bhat <- matrix(0, nrow=J, ncol=R)
  
  non.null.idx <- sort(sample(1:J, 600, replace=FALSE))
  names(non.null.idx) <- names(mu)[non.null.idx]
  
  for(j in 1:J){
    if(j %in% non.null.idx){
      # simulate effect size
      w.j <- sqrt(pmax(abs(rnorm(1, 0, 0.5)), 0.01))
      
      # simulate effect-sharing pattern
      u.j <- ulist[[which(as.numeric(rmultinom(1, 1, pi.u))==1)]]
      beta[j,] <- w.j*ifelse(runif(1) > 0.5, 1, -1)*u.j
    }
    
    bhat[j,] <- rnorm(R, mean=mu[j]+beta[j,], sd=shat[j,])
  }
  
  rownames(beta) <- rownames(shat)
  colnames(beta) <- trts
  rownames(bhat) <- rownames(shat)
  colnames(bhat) <- trts
  
  data <- list(bhat=bhat, shat=shat, beta=beta, mu=mu, non.null.idx=non.null.idx)
  
  saveRDS(data, file=paste0("signal_mash/data/iter", iter, ".Rds"))
}



################################ Plot the true effect-sharing patterns #########################################
library(RColorBrewer)

### define the color map for conditions
cols.trt <- c(brewer.pal(n=9, name="Set1"), brewer.pal(n=8, name="Set2")[-c(2,8)], brewer.pal(n=12, name="Set3")[-c(1,7)])

pdf("effect_sharing_truth.pdf", width=8, height=8)
par(mfrow=c(2,2))

for (k in 1:3){
  v <- ulist[[k]]
  barplot(v/v[which.max(abs(v))], names = trts, cex.names = 0.4, las = 2, col = cols.trt, main = paste0("pattern", k))
}

dev.off()