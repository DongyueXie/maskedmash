
ebci_filename = "spiky0_a6_new_rep"
ash_filename = "ash_Z0_a6_new"
prior = "spiky0"
#ebci_filename = "bimodal0_a6rep"
#ash_filename = "ash_Z0_a6_bimodal"
#prior = "bimodal0"

target = seq(-5,-3,by=0.2)



ebci = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',ebci_filename,'.csv',sep=''))
# number of target is 11, total 5 methods, so for each replicate, there should be 55 counts
table(ebci$id)
# each method should have 11*100 counts
table(ebci$method)

# # first take a look at upper bound for location mixture npmle model
# method = 'amari_locmix_npmle'
#
# par(mfrow=c(3,3))
# for(i in 1:9){
#   idx = (ebci$method==method)&(ebci$id==i)
#   plot(ebci[idx,]$t,ebci[idx,]$upper,type='l')
#   z = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/Z0_a6/Zs_rep',i,'.csv',sep=''))$Column1
#   print(range(z))
# }
#
#
#
# for(i in 1:100){
#   idx = (ebci$method==method)&(ebci$id==i)
#   print(i)
#   print(round(range(ebci[idx,]$lower),3))
# }
mtd_list = c("amari_scalemix","amari_scalemix_eps","amari_scalemix0", "amari_locmix",   "amari_locmix_npmle")
name_list = c('SN-default','SN-eps','SN-0','LN-default','NPMLE')
library(ashr)
if(prior == 'spiky0'){
  g_true = normalmix(c(0.5,0.2,0.1,0.1,0.1),rep(0,5),s=c(0,0.25,0.5,1,2))
}else if(prior == 'bimodal0'){
  g_true = normalmix(c(0.5,0.25,0.25),c(0,-2,2),s=c(0,1,1))
}

bks = seq(0,0.035,length.out = 50)

pdf(file = paste("output/ebci_discussion/",ebci_filename,".pdf",sep=''),
    width = 6,
    height = 4)
# for each method, for each rep, get the z corresponds to lfsr-hat = t, then calc true lfsr of these z-score, plot the histogram.
par(mfrow=c(2,3),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0))
t = 0.01
nreps = 100
for(k in 1:length(mtd_list)){
  zs = c()
  for(i in 1:nreps){
    idx = which((ebci$method==mtd_list[k])&(ebci$id==i))
    if(length(idx)>0){
      temp = ebci[idx,]
      for(j in nrow(temp):2){
        if(temp$upper[j]>t & temp$upper[j-1]<t){
          zs = c(zs,temp$t[j-1]+(temp$t[j]-temp$t[j-1])/(temp$upper[j]-temp$upper[j-1])*(t-temp$upper[j-1]))
          next
        }
      }
    }
  }
  lfsr_true = ash(zs,s=1,g=g_true,fixg=TRUE)$result$lfsr
  hist(lfsr_true,breaks = bks,main=name_list[k],xlab='lfsr',xlim=c(0,3*t),ylim=c(0,50))
  abline(v = t,lty=2)
}
dev.off()


pdf(file = paste("output/ebci_discussion/",ash_filename,"_rep.pdf",sep=''),
    width = 6,
    height = 2)
## ash result
ash_res = readRDS(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',ash_filename,'.rds',sep=''))
ash_lfsr = ash_res$lfsr0
#colnames(ash_lfsr) = ash_res$nullweights
#lfsr_true = ash_res$lfsr_true
par(mfrow=c(1,3),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0))
for(w in 1:length(ash_res$nullweights)){
  zs = c()
  temp = ash_lfsr[,,w]
  for(i in 1:nreps){
      for(j in ncol(temp):2){
        if(temp[i,j]>t & temp[i,j-1]<t){
          zs = c(zs,target[j-1]+(target[j]-target[j-1])/(temp[i,j]-temp[i,j-1])*(t-temp[i,j-1]))
          next
        }
      }
  }
  lfsr_true = ash(zs,s=1,g=g_true,fixg=TRUE)$result$lfsr
  hist(lfsr_true,breaks = bks,main=paste('nullweight =',ash_res$nullweights[w]),xlab='lfsr',xlim=c(0,t*3),ylim = c(0,50))
  abline(v = t,lty=2)
}
dev.off()
