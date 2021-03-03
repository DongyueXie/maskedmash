library(devtools)
load_all('code/mashr/')
files.sources = list.files('code/maskedmashr/')
sapply(files.sources, function(x){source(paste('code/maskedmashr/',x,sep=''))})

library(parallel)
# return a list of : data, and masked mash fit
sim_studyP = function(Ulist,pi,N,seed=12345,nreps = 20,mc.cores = 4,npc=5,df=10){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    # generate data
    simdata = simData.t(N,Ulist,pi,df)
    #datax = mash_set_data(simdata$Bhat,simdata$Shat)
    # fit model
    mash.out = mash_wrapper(simdata$Bhat,P=simdata$P, npc=npc)
    mashonmaskedZ.out = mash_mask(simdata$Bhat,P=simdata$P,npc=npc)
    mashonmaskedZ.out2 = mash_mask(simdata$Bhat,P=simdata$P,thresh=1.5,npc=npc)
    maskedmash.out = maskedmash_wrapper(simdata$Bhat,P=simdata$P,npc=npc)
    maskedmash.out2 = maskedmash_wrapper(simdata$Bhat,P=simdata$P,thresh=1.5,npc=npc)

    list(data = simdata,
         mash.out=mash.out,
         mashonmaskedZ.out=mashonmaskedZ.out,
         mashonmaskedZ.out2=mashonmaskedZ.out2,
         maskedmash.out=maskedmash.out,
         maskedmash.out2=maskedmash.out2)

  },mc.cores = mc.cores)
}



sim_study = function(Ulist,pi,N,seed=12345,nreps = 20,mc.cores = 4,npc=5,df=10){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    # generate data
    simdata = simData.t(N,Ulist,pi,df)
    #datax = mash_set_data(simdata$Bhat,simdata$Shat)
    # fit model
    mash.out = mash_wrapper(simdata$Bhat,npc=npc)
    mashonmaskedZ.out = mash_mask(simdata$Bhat,npc=npc)
    mashonmaskedZ.out2 = mash_mask(simdata$Bhat,thresh=1.5,npc=npc)
    maskedmash.out = maskedmash_wrapper(simdata$Bhat,npc=npc)
    maskedmash.out2 = maskedmash_wrapper(simdata$Bhat,thresh=1.5,npc=npc)

    list(data = simdata,
         mash.out=mash.out,
         mashonmaskedZ.out=mashonmaskedZ.out,
         mashonmaskedZ.out2=mashonmaskedZ.out2,
         maskedmash.out=maskedmash.out,
         maskedmash.out2=maskedmash.out2)

  },mc.cores = mc.cores)
}










#'@title get summary from simu study
#'@description given a sequence of
get.simu.result = function(out,alpha){

  fdr_result = lapply(out,function(x){
    #print(length(x))
    non_null_idx = which(x$data$B!=0)
    #print(non_null_idx)
    fp = lapply(x[-1],function(z){
      rej = mashFDR(x$data$Bhat,z,alpha)$rej.set
      fdps = fdp(rej,non_null_idx)
      powers = powr(rej,non_null_idx)
      c(fdps,powers)
    })

    do.call(cbind,fp)
  })

  temp = do.call(rbind,fdr_result)
  fdr.idx = seq(1,2*length(out),2)
  list(FDP = temp[fdr.idx,],POWER = temp[-fdr.idx,])
}

plot.simu.result = function(obj){
  fdps = obj$FDP
  colnames(fdps) = c("mash",'mash.maskedZ','mash.maskedZ.partial','masked.mash','masked.mash.parital')
  boxplot(fdps,ylab='FDP')
  abline(h = alpha,lty=2)

  powers = obj$POWER
  colnames(powers) = c("mash",'mash.maskedZ','mash.maskedZ.partial','masked.mash','masked.mash.parital')
  boxplot(powers,ylab='Power')
}


######################
######################

signal_sd = sqrt(3)
R = 5
Ulist = list(U1 = matrix(0,nrow=R,ncol=R),
             U2 = diag(R),
             U3 = tcrossprod(c(1,1,0,0,0)),
             U4 = tcrossprod(c(0,0,1,1,1)))
Ulist = lapply(Ulist,function(x){x*signal_sd^2})



out = sim_studyP(Ulist=Ulist,pi = rep(1/length(Ulist),length(Ulist)),
                N = 2000,nreps=30,mc.cores = 4,df=10)
saveRDS(out,file = 'output/maskedmashFDR/t_P_SNR3K4N2000.rds')

####################
####################

out = sim_study(Ulist,pi = rep(1/length(Ulist),length(Ulist)),
                N = 2000,
                nreps=30,mc.cores = 4,df=10)
saveRDS(out,file = 'output/maskedmashFDR/t_SNR3K4N2000.rds')

###################
###################

alpha_list = seq(0,0.5,length.out = 50)
fdr = c()
powers = c()
for(i in 1:length(alpha_list)){
  temp = get.simu.result(out,alpha_list[i])
  fdr = rbind(fdr,colMeans(temp$FDP))
  powers = rbind(powers,colMeans(temp$POWER))
}

#########

plot(alpha_list,fdr[,1],type='l',lwd=2,ylim = range(fdr),xlab = 'Target FDR Level', ylab = 'FDR')
for(i in 2:ncol(fdr)){
  lines(alpha_list,fdr[,i],lwd=2,col=i,lty=i)
}
abline(a=0,b=1,lty=1,col='grey80')
legend('bottomright',c("mash","mash.on.maskedZ.all","mash.on.maskedZ.40%","masked.mash.all","masked.mash.40%"),
       lty = 1:5,col=1:5,lwd=rep(2,5))

#########

plot(alpha_list,powers[,1],type='l',lwd=2,ylim = range(powers),xlab = 'Target FDR Level', ylab = 'Power')
for(i in 2:ncol(powers)){
  lines(alpha_list,powers[,i],lwd=2,col=i,lty=i)
}
legend('bottomright',c("mash","mash.on.maskedZ.all","mash.on.maskedZ.40%","masked.mash.all","masked.mash.40%"),
       lty = 1:5,col=1:5,lwd=rep(2,5))

#########

plot(fdr[,1],powers[,1],type='l',lwd=2,xlim = range(fdr),ylim = range(powers),xlab = 'FDR', ylab = 'Power')
for(i in 2:ncol(powers)){
  lines(fdr[,i],powers[,i],lwd=2,col=i,lty=i)
}
legend('topleft',c("mash","mash.on.maskedZ.all","mash.on.maskedZ.40%","masked.mash.all","masked.mash.40%"),
       lty = 1:5,col=1:5,lwd=rep(2,5))


