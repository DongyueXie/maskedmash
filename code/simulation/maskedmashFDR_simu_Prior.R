library(devtools)
library(parallel)
load_all('code/mashr/')
files.sources = list.files('code/maskedmashr/')
sapply(files.sources, function(x){source(paste('code/maskedmashr/',x,sep=''))})



# return a list of : data, and masked mash fit

#'@param prior normal, t or uniform

sim_study = function(Ulist,pi,N,prior,df,
                     seed=12345,nreps = 20,mc.cores = 4,npc=5,
                     half.uniform=FALSE,mean.range=4){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    # generate data
    simdata = simDataI.ult(N,Ulist,pi,prior,df,half.uniform,mean.range)
    #datax = mash_set_data(simdata$Bhat,simdata$Shat)
    # fit model
    mash.out = mash_wrapper(simdata$Bhat,npc=npc)

    mashonmaskedZ.out = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.5,npc=npc)
    mashonmaskedZ.out2 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.4,npc=npc)
    mashonmaskedZ.out3 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.3,npc=npc)
    mashonmaskedZ.out4 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.2,npc=npc)
    mashonmaskedZ.out5 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.1,npc=npc)
    mashonmaskedZ.out6 = mash_mask(simdata$Bhat,simdata$P,p.thresh = 0.05,npc=npc)

    maskedmash.out = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh = 0.5,npc=npc)
    maskedmash.out2 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.4,npc=npc)
    maskedmash.out3 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.3,npc=npc)
    maskedmash.out4 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.2,npc=npc)
    maskedmash.out5 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.1,npc=npc)
    maskedmash.out6 = maskedmash_wrapper(simdata$Bhat,simdata$P,p.thresh=0.05,npc=npc)

    list(data = simdata,
         mash.out=mash.out,
         mashonmaskedZ.out=mashonmaskedZ.out,
         mashonmaskedZ.out2=mashonmaskedZ.out2,
         mashonmaskedZ.out3=mashonmaskedZ.out3,
         mashonmaskedZ.out4=mashonmaskedZ.out4,
         mashonmaskedZ.out5=mashonmaskedZ.out5,
         mashonmaskedZ.out6=mashonmaskedZ.out6,
         maskedmash.out=maskedmash.out,
         maskedmash.out2=maskedmash.out2,
         maskedmash.out3=maskedmash.out3,
         maskedmash.out4=maskedmash.out4,
         maskedmash.out5=maskedmash.out5,
         maskedmash.out6=maskedmash.out6)

  },mc.cores = mc.cores)

  result

}


######################
######################


R = 5
u = toeplitz(c(1,0.4,0.3,0.2,0.1))
u[4:5,] = 0
u[,4:5] = 0

u2 = toeplitz(c(1,-0.8,0,0,0))
u2[1:3,] = 0
u2[,1:3] = 0

Ulist = list(U1 = matrix(0,nrow=R,ncol=R),
             U2 = diag(R),
             U3 = u,
             U4 = u2)
n = 500



### t prior

# print("running t prior 1")
#
# out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),mean.range = 4,
#                 prior = 't',df = 1,seed=12345,nreps = 30,mc.cores = 4)
# saveRDS(out,file = 'output/maskedmashFDR/t_prior_range4K4N2000df1.rds')
# rm(out)

# print("running t prior 2")
#
# out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),mean.range = 4,
#                 prior = 't',df = 10,seed=12345,nreps = 20,mc.cores = 4)
# saveRDS(out,file = 'output/maskedmashFDR/t_prior_range4K4N2000df10.rds')
# rm(out)
#
# ## unif prior
#
# print("running unif prior 1")
#
# out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),mean.range = 4,
#                 prior = 'uniform',seed=12345,nreps = 20,mc.cores = 4,half.uniform = FALSE)
# saveRDS(out,file = 'output/maskedmashFDR/uniform_prior_range4K4N2000.rds')
# rm(out)
#
# print("running unif prior 2")
#
# out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),mean.range = 2,
#                 prior = 'uniform',seed=12345,nreps = 20,mc.cores = 4,half.uniform = FALSE)
# saveRDS(out,file = 'output/maskedmashFDR/uniform_prior_range2K4N2000.rds')
# rm(out)

## half unif prior

print("running half unif prior 1")

out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),mean.range = 4,
                prior = 'uniform',seed=12345,nreps = 20,mc.cores = 4,half.uniform = TRUE)
saveRDS(out,file = 'output/maskedmashFDR/halfuniform_prior_range4K4N2000.rds')
rm(out)

print("running half unif prior 2")

out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),mean.range = 2,
                prior = 'uniform',seed=12345,nreps = 20,mc.cores = 4,half.uniform = TRUE)
saveRDS(out,file = 'output/maskedmashFDR/halfuniform_prior_range2K4N2000.rds')
rm(out)

### normal prior

print("running normal prior 1")

out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),mean.range = 4,
                prior = 'normal',seed=12345,nreps = 20,mc.cores = 4)
saveRDS(out,file = 'output/maskedmashFDR/normal_prior_range4K4N2000.rds')
rm(out)


print("running normal prior 2")

out = sim_study(Ulist,pi=rep(1/length(Ulist),length(Ulist)),N = n*length(Ulist),mean.range = 3,
                prior = 'normal',seed=12345,nreps = 20,mc.cores = 4)
saveRDS(out,file = 'output/maskedmashFDR/normal_prior_range3K4N2000.rds')
rm(out)

