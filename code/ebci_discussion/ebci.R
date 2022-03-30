
########## generate z scores from spiky #####

spiky = function(n,w = c(0.4,0.2,0.2,0.2),s=c(0.25,0.5,1,2)){
  K = length(w)
  Z = c()
  for (k in 1:K) {
    Z = c(Z,rnorm(n*w[k],0,sqrt(1+s[k]^2)))
  }
  return(Z)
}

##############################################
##############################################

## 1. no point mass when fitting the model

library(ashr)

nreps = 100
n = 5000
target = seq(-5,5,by=0.2)


a = 0.1
mixsd = 0.1
while(a<15.6){
  a = a*1.1
  mixsd = c(mixsd,a)
}

g_true = normalmix(c(0.4,0.2,0.2,0.2),rep(0,4),s=c(0.25,0.5,1,2))

nullweights = c(1,10,100)
PosProb = array(dim=c(nreps,length(target),length(nullweights)))
PosProb_true = matrix(nrow=nreps,ncol=length(target))
set.seed(12345)
for(i in 1:nreps){
  if(i%%5==0){
    print(i)
  }
  z = read.csv(paste('output/ebci/lfsr_no_point_mass/Z/Zs_rep',i,'.csv',sep=''))$Column1
  #z = spiky(n)
  # first get  true lfsr
  PosProb_true[i,] = ash(target,s=1,g=g_true,fixg=TRUE)$result$PositiveProb
  for(j in 1:length(nullweights)){
    ghat = ash(z,s=1,mixcompdist = 'normal',pointmass = FALSE, optmethod='mixSQP',prior = 'nullbiased',mixsd = mixsd,nullweight=nullweights[j])$fitted_g
    res = ash(target,s=1,g=ghat,fixg=TRUE)
    PosProb[i,,j] = res$result$PositiveProb
  }
}

saveRDS(list(PosProb=PosProb,PosProb_true=PosProb_true,g_true=g_true,nullweights=nullweights),file='output/ebci/ash_lfsr_no_point_mass.rds')


## 2. add point mass when fitting the model

library(ashr)

nreps = 100
n = 5000
target = seq(-5,0,by=0.2)
nullweights = c(1,10,100)

a = 0.1
mixsd = 0.1
while(a<15.6){
  a = a*1.1
  mixsd = c(mixsd,a)
}

g_true = normalmix(c(0.4,0.2,0.2,0.2),rep(0,4),s=c(0.25,0.5,1,2))
lfsr = array(dim=c(nreps,length(target),length(nullweights)))
lfsr_true = matrix(nrow=nreps,ncol=length(target))

set.seed(12345)
for(i in 1:nreps){
  if(i%%5==0){
    print(i)
  }
  z = read.csv(paste('output/ebci/lfsr_no_point_mass/Z/Zs_rep',i,'.csv',sep=''))$Column1
  lfsr_true[i,] = ash(target,s=1,g=g_true,fixg=TRUE)$result$lfsr
  for(j in 1:length(nullweights)){
    ghat = ash(z,s=1,mixcompdist = 'normal',pointmass = TRUE, optmethod='mixSQP',
               prior = 'nullbiased',mixsd = mixsd,nullweight=nullweights[j])$fitted_g
    res = ash(target,s=1,g=ghat,fixg=TRUE)
    lfsr[i,,j] = res$result$lfsr
  }
}

saveRDS(list(lfsr=lfsr,lfsr_true=lfsr_true,g_true=g_true,nullweights=nullweights),file='output/ebci/ash_lfsr_pointmass.rds')


## 3. record lfdr


library(ashr)

nreps = 100
n = 5000
target = seq(-5,0,by=0.2)
nullweights = c(1,10,100)

a = 0.1
mixsd = 0.1
while(a<15.6){
  a = a*1.1
  mixsd = c(mixsd,a)
}

g_true = normalmix(c(0.4,0.2,0.2,0.2),rep(0,4),s=c(0.25,0.5,1,2))
lfdr = array(dim=c(nreps,length(target),length(nullweights)))
#lfdr_true = matrix(nrow=nreps,ncol=length(target))

set.seed(12345)
for(i in 1:nreps){
  if(i%%5==0){
    print(i)
  }
  z = read.csv(paste('output/ebci/lfsr_no_point_mass/Z/Zs_rep',i,'.csv',sep=''))$Column1
  #lfdr_true[i,] = ash(target,s=1,g=g_true,fixg=TRUE)$result$lfdr
  for(j in 1:length(nullweights)){
    ghat = ash(z,s=1,mixcompdist = 'normal',pointmass = TRUE, optmethod='mixSQP',
               prior = 'nullbiased',mixsd = mixsd,nullweight=nullweights[j])$fitted_g
    res = ash(target,s=1,g=ghat,fixg=TRUE)
    lfdr[i,,j] = res$result$lfdr
  }
}

saveRDS(list(lfdr=lfdr,g_true=g_true,nullweights=nullweights),file='output/ebci/ash_lfdr.rds')

















plot(target,colMeans(PosProb[,,1]),type='l',xlab = 'z',ylab = expression(paste('P(',mu>=0,'|',Z==z,')')))
for(j in 2:length(nullweights)){
  lines(target,colMeans(PosProb[,,j]),col=j)
}
legend('bottomright',paste('nullweight=',nullweights,sep=''),lty=rep(1,5),col=1:5)

# Zoom in
idx = 1:16
plot(target[idx],colMeans(PosProb[,idx,1]),type='l',xlab = 'z',ylab = expression(paste('P(',mu>=0,'|',Z==z,')')))
for(j in 2:length(nullweights)){
  lines(target[idx],colMeans(PosProb[,idx,j]),col=j)
}
legend('topleft',paste('nullweight=',nullweights,sep=''),lty=rep(1,5),col=1:5)

# Zoom in
idx = 1:11
plot(target[idx],colMeans(PosProb[,idx,1]),type='l',xlab = 'z',ylab = expression(paste('P(',mu>=0,'|',Z==z,')')))
for(j in 2:length(nullweights)){
  lines(target[idx],colMeans(PosProb[,idx,j]),col=j)
}
abline(h = 0.01,lty=2)
legend('topleft',paste('nullweight=',nullweights,sep=''),lty=rep(1,5),col=1:5)

# Zoom in
idx = 1:6
plot(target[idx],colMeans(PosProb[,idx,1]),type='l',xlab = 'z',ylab = expression(paste('P(',mu>=0,'|',Z==z,')')))
for(j in 2:length(nullweights)){
  lines(target[idx],colMeans(PosProb[,idx,j]),col=j)
}
legend('topleft',paste('nullweight=',nullweights,sep=''),lty=rep(1,5),col=1:5)


# add point mass
lfsr = array(dim=c(nreps,length(target),length(nullweights)))
set.seed(12345)
for(i in 1:nreps){
  if(i%%5==0){
    print(i)
  }
  z = spiky(n)
  for(j in 1:length(nullweights)){
    ghat = ash(z,s=1,mixcompdist = 'normal',pointmass = TRUE, optmethod='mixSQP',prior = 'nullbiased',mixsd = mixsd,nullweight=nullweights[j])$fitted_g
    res = ash(target,s=1,g=ghat,fixg=TRUE)
    lfsr[i,,j] = res$result$lfsr
  }
}

idx = 1:26
plot(target[idx],colMeans(lfsr[,idx,1]),type='l',xlab = 'z',ylab = expression(paste('P(',mu>=0,'|',Z==z,')')),ylim=range(lfsr))
for(j in 2:length(nullweights)){
  lines(target[idx],colMeans(lfsr[,idx,j]),col=j)
}
legend('topleft',paste('nullweight=',nullweights,sep=''),lty=rep(1,5),col=1:5)


idx = 1:16
plot(target[idx],colMeans(lfsr[,idx,1]),type='l',xlab = 'z',ylab = expression(paste('P(',mu>=0,'|',Z==z,')')),ylim=range(lfsr[,idx,]))
for(j in 2:length(nullweights)){
  lines(target[idx],colMeans(lfsr[,idx,j]),col=j)
}
legend('topleft',paste('nullweight=',nullweights,sep=''),lty=rep(1,5),col=1:5)

idx = 1:11
plot(target[idx],colMeans(lfsr[,idx,1]),type='l',xlab = 'z',ylab = expression(paste('P(',mu>=0,'|',Z==z,')')),ylim=range(lfsr[,idx,]))
for(j in 2:length(nullweights)){
  lines(target[idx],colMeans(lfsr[,idx,j]),col=j)
}
abline(h=0.01,lty=2)
legend('topleft',paste('nullweight=',nullweights,sep=''),lty=rep(1,5),col=1:5)

idx = 1:6
plot(target[idx],colMeans(lfsr[,idx,1]),type='l',xlab = 'z',ylab = expression(paste('P(',mu>=0,'|',Z==z,')')),ylim=range(lfsr[,idx,]))
for(j in 2:length(nullweights)){
  lines(target[idx],colMeans(lfsr[,idx,j]),col=j)
}
legend('topleft',paste('nullweight=',nullweights,sep=''),lty=rep(1,5),col=1:5)
