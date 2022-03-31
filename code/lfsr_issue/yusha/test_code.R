library(MASS)
library(devtools)
load_all('code/mashr/')
###################################
######### generate data ###########
###################################

#'@description generate bhat given shat
#'@param shat
#'@return a list of bhat,shat,b
generate_data = function(Ulist,pi.u,shat,mu,non_null_n = 600){
  library(MASS)

  gene_names = rownames(shat)
  cond_names = colnames(shat)

  n = nrow(shat)
  p = ncol(shat)
  K = length(Ulist)

  b = matrix(nrow=n,ncol = p)
  bhat = matrix(nrow = n,ncol=p)

  non_null_idx = sample(1:n,non_null_n)
  non_null_idx_keep = non_null_idx
  #generate non-null samples
  for(k in 1:K){
    w = pi.u[k]
    nk = round(non_null_n*w)
    idxk = non_null_idx[1:nk]
    non_null_idx = non_null_idx[-(1:nk)]

    bk = matrix(nrow=nk,ncol=p)
    bhatk = matrix(nrow=nk,ncol=p)
    for(i in 1:nk){
      bk[i,] = mvrnorm(1,rep(mu[idxk[i]],p),tcrossprod(Ulist[[k]]))
      bhatk[i,] = rnorm(p,mean = bk[i,],sd = shat[idxk[i],])
    }
    b[idxk,] = bk
    bhat[idxk,] = bhatk
  }

  # generate null samples
  b_null = mu[-non_null_idx_keep]%*%t(rep(1,p))
  bhat_null = matrix(rnorm((n-non_null_n)*p,b_null,shat[-non_null_idx_keep,]),nrow = n-non_null_n,ncol=p)

  b[-non_null_idx_keep,] = b_null
  bhat[-non_null_idx_keep,] = bhat_null

  rownames(b) = gene_names
  rownames(bhat) = gene_names

  colnames(b) = cond_names
  colnames(bhat) = cond_names

  return(list(b=b,bhat=bhat,shat=shat,
              Ulist=Ulist,pi.u=pi.u,
              shat=shat,mu=mu,non_null_n=non_null_n))
}


sim_design <- readRDS("code/lfsr_issue/yusha/sim_design.rds")
library(mashr)


######################################################
################ check the calibration ###############
######################################################

# calculate nominal FSR and empirical FSR, as well as empirical FDR.

fdp = function(rej.idx, null.idx){
  if(length(rej.idx)==0){
    return(0)
  }else{
    sum(rej.idx%in%null.idx)/length(rej.idx)
  }
}

fsp = function(rej.idx,est.sign,true.sign){
  if(length(rej.idx)==0){
    return(0)
  }else{
    sum(est.sign[rej.idx]!=true.sign[rej.idx])/length(rej.idx)
  }
}

#'@title mash for FDR control
#'@description order the tests from most to least significant based on lfsr, and reject least significant ones until FDP<alpha
#'@param lfsr vectorized lfsr hat from mash
#'@param alpha target FSR level

mashFSR = function(lfsr,alpha = 0.05){

  # rank the tests from most significant to least significant
  # and also as initial rejection set
  if(mean(lfsr)<=alpha){
    rej.set = 1:length(lfsr)
  }else{
    l = order(lfsr,decreasing = FALSE)
    lfsr.ordered  =lfsr[l]
    lfsr.cmean = cumsum(lfsr.ordered)/(1:length(lfsr.ordered))
    rej.set = l[which(lfsr.cmean<=alpha)]
  }

  return(rej.set)

}


######################################################
######### let's first study the most basic case ###########
######################################################
# 1. all mu are 0
# 2. run 20 replicates
# 3. look at the nominal FSR and empirical FSR,and FDR.
# 4. for each replicate, also look at lfsr_hat vs true_lfsr



Ulist = sim_design$ulist
pi.u = sim_design$pi.u
shat = sim_design$shat
#mu = sim_design$mu
mu = rep(0,length(sim_design$mu))

set.seed(12345)
# use a smaller set for this study
ii = sample(1:nrow(shat),2000)
shat = shat[ii,]
mu = mu[ii]

res = list()

set.seed(12345)
for(idx in 1:10){
  ################################ generate data #########################################
  print(sprintf("Running replicate %d", idx))
  print('Generating simulation dataset')

  datax = generate_data(Ulist=Ulist,pi.u = pi.u,shat=shat,mu = mu)
  #data <- readRDS(file=paste0("data/iter", idx, ".Rds"))


  ################################ Run mash by directly subtracting the median ##################################
  start_time = proc.time()

  b = datax$b
  bhat <- datax$bhat
  shat <- datax$shat
  R <- ncol(bhat)
  b = b - datax$mu %*% t(rep(1,R))
  bhat <- bhat - datax$mu %*% t(rep(1,R))

  ### set up for mash
  mash.data = mash_set_data(bhat, shat, alpha = 1)
  ################################################################################################
  ################################################################################################

  ### true prior covariances
  U = list()
  U[[1]] <- Ulist[[1]] %*% t(Ulist[[1]])
  U[[2]] <- Ulist[[2]] %*% t(Ulist[[2]])
  U[[3]] <- Ulist[[3]] %*% t(Ulist[[3]])
  names(U) <- c("PC1", "PC2", "PC3")

  ### fit mash

  # ### use true cov
  #
  # mash.fit = mash(mash.data, U, algorithm.version = 'Rcpp',verbose=F)
  # runtime = proc.time() - start_time
  # mash.fit[["runtime"]] = runtime

  ### use ed estimated ones
  print('running ED')
  U.c = cov_canonical(mash.data)
  m.1by1 = mash_1by1(mash.data, alpha=1)
  strong = get_significant_results(m.1by1, 0.2)
  U.pca = cov_pca(mash.data, 5, subset=strong)
  U.ed = cov_ed(mash.data, U.pca, subset=strong)
  print('fitting mash')
  mash.fit = mash(mash.data, c(U.c,U.ed), algorithm.version = 'Rcpp',verbose=F)

  res[[idx]] = list(mash.fit = mash.fit,data = datax)

}


####################################
############ get result ############
####################################

get_summary = function(res,target_fsr){

  empFSR = c()
  empFDR = c()

  for(i in 1:length(target_fsr)){
    out = lapply(res,function(z){
      rej.idx = mashFSR(c(z$mash.fit$result$lfsr),target_fsr[i])
      null.idx = which(z$data$b==0)
      est.sign = sign(c(z$mash.fit$result$PosteriorMean))
      true.sign = sign(c(z$data$b))
      return(c(fsp(rej.idx,est.sign,true.sign),fdp(rej.idx,null.idx)))
    })
    out = do.call(rbind,out)
    #print(out)
    empFSR = c(empFSR,colMeans(out)[1])
    empFDR = c(empFDR,colMeans(out)[2])
  }

  return(list(empFSR=empFSR,empFDR=empFDR))

}

target_fsr = c(0.01,0.03,0.05,0.07,0.1,0.15,0.2,0.25)
out = get_summary(res,target_fsr)
plot(target_fsr,out$empFSR,type = 'l',ylab='empirical FSR',xlab='target FSR')
abline(a=0,b=1,lty=2)
plot(target_fsr,out$empFDR,type = 'l',ylab='empirical FDR',xlab='target FSR')
abline(a=0,b=1,lty=2)







