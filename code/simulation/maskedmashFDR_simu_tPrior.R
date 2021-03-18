library(devtools)
library(parallel)
load_all('code/mashr/')
files.sources = list.files('code/maskedmashr/')
sapply(files.sources, function(x){source(paste('code/maskedmashr/',x,sep=''))})



# return a list of : data, and masked mash fit
sim_study = function(Ulist,pi,N,prior,df,signal_sd = sqrt(3),
                     seed=12345,nreps = 20,mc.cores = 4,npc=5,
                     half.uniform=FALSE,unif.range=3){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    # generate data
    simdata = simDataI.ult(N,Ulist,pi,prior,signal_sd,df,half.uniform,unif.range)
    #datax = mash_set_data(simdata$Bhat,simdata$Shat)
    # fit model
    mash.out = mash_wrapper(simdata$Bhat,npc=npc)
    mashonmaskedZ.out = mash_mask(simdata$Bhat,npc=npc)
    #mashonmaskedZ.out2 = mash_mask(simdata$Bhat,thresh=thresh,npc=npc)
    maskedmash.out = maskedmash_wrapper(simdata$Bhat,npc=npc)
    maskedmash.out2 = maskedmash_wrapper(simdata$Bhat,thresh=1,npc=npc)
    maskedmash.out3 = maskedmash_wrapper(simdata$Bhat,thresh=1.5,npc=npc)
    maskedmash.out4 = maskedmash_wrapper(simdata$Bhat,thresh=2,npc=npc)

    list(data = simdata,
         mash.out=mash.out,
         mashonmaskedZ.out=mashonmaskedZ.out,
         #mashonmaskedZ.out2=mashonmaskedZ.out2,
         maskedmash.out=maskedmash.out,
         maskedmash.out2=maskedmash.out2,
         maskedmash.out3=maskedmash.out3,
         maskedmash.out4=maskedmash.out4)

  },mc.cores = mc.cores)

  result

}

#' #'@title get summary from simu study
#' #'@description given a sequence of
#' get.simu.result = function(out,alpha){
#'
#'   fdr_result = mclapply(out,function(x){
#'     #print(length(x))
#'     non_null_idx = which(x$data$B!=0)
#'     #print(non_null_idx)
#'     rej.mash = mashFDR(x[[2]],alpha)
#'     fp.mash = c(fdp(rej.mash,non_null_idx),powr(rej.mash,non_null_idx))
#'     fp = lapply(x[-c(1:2)],function(z){
#'       rej = maskedmashFDR(x$data$Bhat,z,alpha)$rej.set
#'       fdps = fdp(rej,non_null_idx)
#'       powers = powr(rej,non_null_idx)
#'       c(fdps,powers)
#'     })
#'
#'     cbind(fp.mash,do.call(cbind,fp))
#'   },mc.cores = 4)
#'
#'   temp = do.call(rbind,fdr_result)
#'   fdr.idx = seq(1,2*length(out),2)
#'   list(FDP = temp[fdr.idx,],POWER = temp[-fdr.idx,])
#' }


######################
######################


R = 5
u = toeplitz(c(1,0.4,0.3,0.2,0.1))
u[3:5,] = 0
u[,3:5] = 0

u2 = toeplitz(c(1,-0.4,0.3,0.2,0.1))
u2[1:3,] = 0
u2[,1:3] = 0

Ulist = list(U1 = matrix(0,nrow=R,ncol=R),
             U2 = diag(R),
             U3 = u,
             U4 = u2)
n = 500

signal_sd = sqrt(2)
out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),
                prior = 't',df = 5,signal_sd = signal_sd,seed=12345,nreps = 20,mc.cores = 4,npc=3)
saveRDS(out,file = 'output/maskedmashFDR/t_prior_SNR2K4N2000df5.rds')
rm(out)


out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),prior = 'uniform',
                seed=12345,nreps = 20,mc.cores = 4,npc=3,half.uniform = FALSE,unif.range = 3)
saveRDS(out,file = 'output/maskedmashFDR/uniform_prior_range3K4N2000.rds')

out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),prior = 'uniform',
                seed=12345,nreps = 20,mc.cores = 4,npc=3,half.uniform = TRUE,unif.range = 3)
saveRDS(out,file = 'output/maskedmashFDR/halfuniform_prior_range3K4N2000.rds')
