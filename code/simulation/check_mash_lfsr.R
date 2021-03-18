library(devtools)
library(parallel)
load_all('code/mashr/')
files.sources = list.files('code/maskedmashr/')
sapply(files.sources, function(x){source(paste('code/maskedmashr/',x,sep=''))})



## look at mash calibration of lfsr

simu_mash_lfsr = function(Ulist,pi,N,signal_sd = sqrt(3),seed=12345,
                          nreps = 20,mc.cores = 4,npc=5,
                          use.U.true = FALSE,adjust='lb',
                          prior='normal',df=10,half.uniform=FALSE,unif.range=3){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    # generate data
    simdata = simDataI.ult(N,Ulist,pi,prior=prior,signal_sd,
                           df=df,half.uniform=half.uniform,unif.range=unif.range)
    #datax = mash_set_data(simdata$Bhat,simdata$Shat)
    # fit model
    if(use.U.true){
      mash.out = mash_wrapper(simdata$Bhat,npc=npc,U.true = Ulist)
    }else{
      mash.out = mash_wrapper(simdata$Bhat,npc=npc,adjust=adjust)
    }


    list(data = simdata,
         mash.out=mash.out)

  },mc.cores = mc.cores)

  result

}



get.fsr = function(out,alpha){

  fdr_result = mclapply(out,function(x){
    #print(length(x))
    non_null_idx = which(x$data$B!=0)
    #print(non_null_idx)
    rej.mash = mashFDR(x$mash.out,alpha)
    fp.mash = c(fdp(rej.mash,non_null_idx),fsp(x$data$Bhat[rej.mash],x$data$B[rej.mash]),powr(rej.mash,non_null_idx))

    fp.mash
  },mc.cores = 4)

  temp = do.call(rbind,fdr_result)

  colnames(temp) = c('FDR','FSR','power')

  temp
}

plot.out = function(out,
                    alpha_list= seq(0,0.3,length.out = 30),
                    titles = NULL){
  fsrs = c()
  fdrs = c()
  powers = c()
  for(i in 1:length(alpha_list)){
    temp = get.fsr(out,alpha_list[i])
    fsrs[i] = colMeans(temp)[2]
    fdrs[i] = colMeans(temp)[1]
    powers[i] = colMeans(temp)[3]
  }
  plot(alpha_list,fsrs,type='l',lwd=2,ylim = range(fsrs),xlab = 'Target FSR Level', ylab = 'FSR',main = titles)
  abline(a=0,b=1,lty=1,col='grey80')

  plot(alpha_list,fdrs,type='l',lwd=2,ylim = range(fdrs),xlab = 'Target FDR Level', ylab = 'FDR',main = titles)
  abline(a=0,b=1,lty=1,col='grey80')

  plot(alpha_list,powers,type='l',lwd=2,ylim = range(powers),xlab = 'Target FDR Level', ylab = 'POWER',main = titles)
}


# plot.fdr = function(out,
#                     alpha_list= seq(0,0.3,length.out = 30),
#                     titles = NULL){
#   #fsrs = c()
#   fdrs = c()
#   #powers = c()
#   for(i in 1:length(alpha_list)){
#     temp = get.fsr(out,alpha_list[i])
#     #fsrs[i] = colMeans(temp)[2]
#     fdrs[i] = colMeans(temp)[1]
#     #powers[i] = colMeans(temp)[3]
#   }
#   plot(alpha_list,fdrs,type='l',lwd=2,ylim = range(fdrs),xlab = 'Target FDR Level', ylab = 'FDR',main = titles)
#   abline(a=0,b=1,lty=1,col='grey80')
# }

