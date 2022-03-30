
# run ash, spiky prior

prior = 'spiky'
z_folder = 'Z_a6'
ash_filename = 'ash_Z_a6_new'

library(ashr)
target = seq(-5,-3,by=0.2)
nullweights = c(1,10,100)
nreps = 100
n = 5000

a = 0.1
mixsd = 0.1
while(a<15.6){
  a = a*1.1
  mixsd = c(mixsd,a)
}

if(prior == 'spiky'){
  g_true = normalmix(2*c(0.2,0.1,0.1,0.1),rep(0,4),s=c(0.25,0.5,1,2))
}

if(prior == 'spiky0'){
  g_true = normalmix(c(0.5,0.2,0.1,0.1,0.1),rep(0,5),s=c(0,0.25,0.5,1,2))
}

lfsr = array(dim=c(nreps,length(target),length(nullweights)))
lfsr0 = array(dim=c(nreps,length(target),length(nullweights)))
lfsr_true = ash(target,s=1,g=g_true,fixg=TRUE)$result$lfsr

set.seed(12345)
for(i in 1:nreps){
  if(i%%5==0){
    print(i)
  }
  z = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',z_folder,'/Zs_rep',i,'.csv',sep=''))$Column1
  for(j in 1:length(nullweights)){
    ghat = ash(z,s=1,mixcompdist = 'normal',pointmass = TRUE, optmethod='mixSQP',
               prior = 'nullbiased',mixsd = mixsd,nullweight=nullweights[j])$fitted_g
    res = ash(target,s=1,g=ghat,fixg=TRUE)
    lfsr0[i,,j] = res$result$lfsr

  }
}

saveRDS(list(lfsr0=lfsr0,lfsr_true=lfsr_true,
             g_true=g_true,nullweights=nullweights),
        file=paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',ash_filename,'.rds',sep=''))



# run ash, bimodal prior

library(ashr)
target = seq(-5,-3,by=0.2)
nullweights = c(1,10,100)
nreps = 100
n = 5000

a = 0.1
mixsd = 0.1
while(a<15.6){
  a = a*1.1
  mixsd = c(mixsd,a)
}

#g_true = normalmix(c(0.5,0.25,0.25),c(0,-2,2),s=c(0,1,1))
g_true = normalmix(c(0.5,0.5),c(-2,2),s=c(1,1))
#g_true = normalmix(2*c(0.2,0.1,0.1,0.1),rep(0,4),s=c(0.25,0.5,1,2))
#lfsr = array(dim=c(nreps,length(target),length(nullweights)))
lfsr0 = array(dim=c(nreps,length(target),length(nullweights)))
lfsr_true = ash(target,s=1,g=g_true,fixg=TRUE)$result$lfsr

set.seed(12345)
for(i in 1:nreps){
  if(i%%5==0){
    print(i)
  }
  z = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/Z_a6_bimodal/Zs_rep',i,'.csv',sep=''))$Column1
  for(j in 1:length(nullweights)){
    ghat = ash(z,s=1,mixcompdist = 'normal',pointmass = TRUE, optmethod='mixSQP',
               prior = 'nullbiased',mixsd = mixsd,nullweight=nullweights[j])$fitted_g
    res = ash(target,s=1,g=ghat,fixg=TRUE)
    lfsr0[i,,j] = res$result$lfsr
  }
}

saveRDS(list(lfsr0=lfsr0,lfsr_true=lfsr_true,
             g_true=g_true,nullweights=nullweights),
        file='/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/ash_Z_a6_bimodal.rds')
