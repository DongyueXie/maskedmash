# check each reps

ebci_filename = "spiky_a6_new_rep"
ebci = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',ebci_filename,'.csv',sep=''))

method = "amari_locmix_npmle"
par(mfrow=c(4,5))
for(i in 1:100){
  idx = which((ebci$method==method)&(ebci$id==i))
  if(length(idx)>0){
    plot(ebci[idx,]$t,ebci[idx,]$lower,type='l',main=i,ylim=c(0,0.05))
    z = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/Z_a6/Zs_rep',i,'.csv',sep=''))$Column1
    print(range(z))
  }
}

i = 59
idx = which((ebci$method==method)&(ebci$id==i))
par(mfrow=c(1,1))
plot(ebci[idx,]$t,ebci[idx,]$lower,type='l',main=i,ylim=c(0,10))

z = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/Z_a6/Zs_rep',i,'.csv',sep=''))$Column1
print(range(z))
hist(z,breaks = 100)


# we manually remove the 59th interations result from spiky NPMLE
ebci_save = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/','spiky_a6_new','.csv',sep=''))

ebci_save[ebci_save$method=='amari_locmix_npmle',]$lower_mean[6] = mean(ebci[(ebci$method=='amari_locmix_npmle')&(ebci$t==-4)&(ebci$id!=59),]$lower)

write.csv(ebci_save,file='/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/spiky_a6_new_rm59.csv')
