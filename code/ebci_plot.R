
ebci = read.csv('output/ebci_discussion/spiky.csv')
ash_res = readRDS('output/ebci_discussion/ash_Z_a6.rds')
ash_lfsr = apply(ash_res$lfsr0,3,colMeans)
colnames(ash_lfsr) = ash_res$nullweights
lfsr_true = ash_res$lfsr_true


pdf(file = "output/ebci_discussion/test25.pdf",
    width = 16,
    height = 6)

par(mfrow=c(2,5))
for(method in mtd_list){
  ploter(method,ebci,ash_lfsr,lfsr_true,idx=1:11,target=seq(-5,-3,by=0.2))
}


ash_res = readRDS('output/ebci_discussion/ash_Z0_a6.rds')
ash_lfsr = apply(ash_res$lfsr0,3,colMeans)
colnames(ash_lfsr) = ash_res$nullweights
lfsr_true = ash_res$lfsr_true
ebci = read.csv('output/ebci_discussion/spiky0_a6.csv')

for(method in mtd_list){
  ploter(method,ebci,ash_lfsr,lfsr_true,idx=1:11,target=seq(-5,-3,by=0.2),ylim=c(0,0.5))
}

dev.off()

