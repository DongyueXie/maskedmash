ebci_filename0 = "spiky0_a6_new"
ebci_filename = "spiky_a6_new_rm59"
ash_filename = 'ash_Z_a6_new'
ash_filename0 = 'ash_Z0_a6_new'
outname = 'spiky_a6_new_revise3'

ploter = function(method,ebci,ash_lfsr,lfsr_true,idx,target,color = rgb(0, 0, 0.8,0.4),ylim=c(0,0.5),nm,ploti,main = '',bty='l'){
  out = ebci[ebci$method == method,]
  if(ploti ==5){
    plot(target[idx],lfsr_true[idx],type='l',xlab='',
         ylab=expression(paste('P(',mu>=0,'|',Z==z,')')),
         ylim=ylim,
         bty=bty)
  }else if(ploti==1){
    plot(target[idx],lfsr_true[idx],type='l',xlab='',
         ylab=expression(paste('P(',mu>=0,'|',Z==z,')')),
         ylim=ylim,
         main=main,
         bty=bty)
  }else{
    plot(target[idx],lfsr_true[idx],type='l',xlab='',
         ylab=expression(paste('P(',mu>=0,'|',Z==z,')')),
         ylim=ylim,
         bty=bty)
  }

  #axis(2, at = c(0,0.05,0.1,0.2,0.3,0.4),labels = c(0,0.05,0.1,0.2,0.3,0.4))
  polygon(c(rev(out$t[idx]), out$t[idx]), c(rev(out$upper_mean[idx]), out$lower_mean[idx]), col=color, border = FALSE)
  lines(target[idx],lfsr_true[idx],col=2)
  ltys=c(5,3,4)
  for(j in 1:3){
    lines(target[idx],ash_lfsr[idx,j],lty=ltys[j])
  }
  abline(h=0.01,lty=2,col='grey80')
  abline(h=0.05,lty=2,col='grey80')
  if(ploti==1){
    legend('topleft',c('Ground truth',paste('nullweight = ',ash_res$nullweights),'AMARI CI'),
           bty="n",
           lty=c(1,5,3,4,0),
           lwd=c(rep(1,4),0),
           col = c(2,1,1,1, color),
           #fill=c(rep(NA,4), color),
           #border = rep(NA,5),
           pch = c(rep(NA,4),22),
           pt.bg = c(rep(NA,4),color),
           pt.cex = 2)
  }
}
mtd_list = c("amari_scalemix","amari_scalemix_eps","amari_scalemix0", "amari_locmix",   "amari_locmix_npmle")
name_list = c('SN default','SN eps','SN 0','LN default','LN NPMLE')

ebci = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',ebci_filename,'.csv',sep=''))
ash_res = readRDS(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',ash_filename,'.rds',sep=''))
ash_lfsr = apply(ash_res$lfsr0,3,colMeans)
colnames(ash_lfsr) = ash_res$nullweights
lfsr_true = ash_res$lfsr_true


#pdf(file = paste("output/ebci_discussion/",outname,".pdf",sep=''),
#    width = 6,
#    height = 12)

setEPS()
postscript("output/ebci_discussion/try.eps", width = 6,height = 12)

par(mfcol=c(5,2),mar=c(1.5, 2.5, 0.2, 0.4),mgp=c(1.4, 0.5, 0))
for(i in 1:length(mtd_list)){
  ploter(mtd_list[i],ebci,ash_lfsr,lfsr_true,idx=1:11,target=seq(-5,-3,by=0.2),ylim=c(0,0.4),nm=name_list[i],ploti = i)
}


ploter_right = function(method,ebci,ash_lfsr,lfsr_true,idx,target,color = rgb(0, 0, 0.8,0.4),ylim=c(0,0.5),nm,ploti,main='',bty='l'){
  out = ebci[ebci$method == method,]
  if(ploti==5){
    plot(target[idx],lfsr_true[idx],type='l',xlab='',
         ylab='',
         ylim=ylim,
         bty=bty)
  }else if(ploti==1){
    plot(target[idx],lfsr_true[idx],type='l',xlab='',
         ylab='',
         ylim=ylim,
         main=main,
         bty=bty)
  }else{
    plot(target[idx],lfsr_true[idx],type='l',xlab='',
         ylab='',
         ylim=ylim,
         bty=bty)
  }
  #axis(2, at = c(0,0.05,0.1,0.2,0.3,0.4),labels = c(0,0.05,0.1,0.2,0.3,0.4))
  polygon(c(rev(out$t[idx]), out$t[idx]), c(rev(out$upper_mean[idx]), out$lower_mean[idx]), col=color, border = FALSE)
  lines(target[idx],lfsr_true[idx],col=2)
  ltys=c(5,3,4)
  for(j in 1:3){
    lines(target[idx],ash_lfsr[idx,j],lty=ltys[j])
  }
  abline(h=0.01,lty=2,col='grey80')
  abline(h=0.05,lty=2,col='grey80')
  # legend('topleft',c('Ground truth',paste('nullweight = ',ash_res$nullweights),nm),
  #        bty="n",
  #        lty=c(1,5,3,4,0),
  #        lwd=c(rep(1,4),0),
  #        col = c(2,1,1,1, color),
  #        #fill=c(rep(NA,4), color),
  #        #border = rep(NA,5),
  #        pch = c(rep(NA,4),22),
  #        pt.bg = c(rep(NA,4),color),
  #        pt.cex = 2)
}

ebci = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',ebci_filename0,'.csv',sep=''))
ash_res = readRDS(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',ash_filename0,'.rds',sep=''))
ash_lfsr = apply(ash_res$lfsr0,3,colMeans)
colnames(ash_lfsr) = ash_res$nullweights
lfsr_true = ash_res$lfsr_true


for(i in 1:length(mtd_list)){
  ploter_right(mtd_list[i],ebci,ash_lfsr,lfsr_true,idx=1:11,target=seq(-5,-3,by=0.2),ylim=c(0,0.4),nm=name_list[i],ploti=i)
}

dev.off()

# ##############################
# ######## use ggplot ##########
# ##############################
# ebci_filename0 = "spiky0_a6_new"
# ebci_filename = "spiky_a6_new_rm59"
# ash_filename = 'ash_Z_a6_new'
# ash_filename0 = 'ash_Z0_a6_new'
# outname = 'spiky_a6_new_ggplot'
#
# mtd_list = c("amari_scalemix","amari_scalemix_eps",
#              "amari_scalemix0", "amari_locmix",
#              "amari_locmix_npmle")
# name_list = c('SN default','SN eps','SN 0','LN default','NPMLE')
# library(reshape2)
# library(ggplot2)
# library(gridExtra)
# ploter_ggplot = function(method,ebci,ash_lfsr,lfsr_true,target = seq(-5,-3,by=0.2),nm){
#   rownames(ash_lfsr) = target
#   ash_lfsr = melt(ash_lfsr,varnames=c('z','nullweight'),value.name = 'lfsr')
#   ash_lfsr$nullweight = as.factor(ash_lfsr$nullweight)
#   lfsr_true=data.frame(z=target,lfsr=lfsr_true)
#
#   pp = ggplot(ebci[ebci$method==method,],aes(x=t))+
#     geom_ribbon(aes(ymin = lower_mean, ymax = upper_mean, fill = nm)) +ylim(0, 0.4)+scale_fill_manual("",values=rgb(0, 0, 0.8,0.4))+
#     geom_line(data=lfsr_true,aes(x=z,y=lfsr,color='Ground truth'))+
#     geom_line(data=ash_lfsr,aes(x=z,y=lfsr,group=nullweight,linetype=nullweight))+
#     scale_linetype_manual(values=c("longdash", "dotted", "dotdash"))+
#     scale_color_manual("",values=c("red"))+
#     geom_hline(yintercept = 0.05,linetype="dashed",color="red",alpha=0.5)+
#     labs(x = "z")
#   pp
# }
#
#
# ebci = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',
#                       ebci_filename,'.csv',sep=''))
# ebci$lower_mean = pmax(ebci$lower_mean,0)
# ash_res = readRDS(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',
#                         ash_filename,'.rds',sep=''))
# ash_lfsr = apply(ash_res$lfsr0,3,colMeans)
# colnames(ash_lfsr) = ash_res$nullweights
# lfsr_true = ash_res$lfsr_true
#
#
#
#
# p_list1 = list()
# p_list1[[1]]=ploter_ggplot(mtd_list[1],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[1])
# p_list1[[2]]=ploter_ggplot(mtd_list[2],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[2])
# p_list1[[3]]=ploter_ggplot(mtd_list[3],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[3])
# p_list1[[4]]=ploter_ggplot(mtd_list[4],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[4])
# p_list1[[5]]=ploter_ggplot(mtd_list[5],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[5])
#
#
#
#
# ebci = read.csv(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',
#                       ebci_filename0,'.csv',sep=''))
# ash_res = readRDS(paste('/project2/mstephens/dongyue/empirical-bayes-confidence-intervals-paper/',
#                         ash_filename0,'.rds',sep=''))
# ash_lfsr = apply(ash_res$lfsr0,3,colMeans)
# colnames(ash_lfsr) = ash_res$nullweights
# lfsr_true = ash_res$lfsr_true
#
#
# p_list2 = list()
# p_list2[[1]]=ploter_ggplot(mtd_list[1],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[1])
# p_list2[[2]]=ploter_ggplot(mtd_list[2],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[2])
# p_list2[[3]]=ploter_ggplot(mtd_list[3],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[3])
# p_list2[[4]]=ploter_ggplot(mtd_list[4],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[4])
# p_list2[[5]]=ploter_ggplot(mtd_list[5],ebci,ash_lfsr,lfsr_true,target=seq(-5,-3,by=0.2),nm=name_list[5])
#
#
# Plotlist = c(p_list1,p_list2)
# glist = lapply(Plotlist, ggplotGrob)
# ggsave(
#   filename = paste("output/ebci_discussion/",outname,".pdf",sep=''),
#   plot = marrangeGrob(glist, nrow=5, ncol=2,top=NULL),
#   width = 7, height = 12
# )
#
#
