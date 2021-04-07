library(devtools)
library(parallel)
library(mvtnorm)
library(udr)
load_all('code/mashr/')
files.sources = list.files('code/maskedmashr/')
sapply(files.sources, function(x){source(paste('code/maskedmashr/',x,sep=''))})



## look at mash calibration of lfsr

simu_mash_lfsr = function(Ulist,pi,N,seed=12345,mean.range = 4,
                          nreps = 20,mc.cores = 4,npc=5,
                          use.U.true = FALSE,adjust='lb',
                          prior='normal',df=10,half.uniform=FALSE){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    # generate data
    simdata = simDataI.ult(N,Ulist,pi,prior=prior,
                           df=df,half.uniform=half.uniform,mean.range=mean.range)
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


simu_mash_lfsr_search = function(Ulist,pi,N,s,
                                 seed=12345,mean.range = 4,
                          nreps = 20,mc.cores = 4,npc=5,
                          prior='normal',df=10,half.uniform=FALSE){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    # generate data
    simdata = simDataI.ult(N,Ulist,pi,prior=prior,
                           df=df,half.uniform=half.uniform,mean.range=mean.range)
    #datax = mash_set_data(simdata$Bhat,simdata$Shat)
    # fit model
    mash.out = mash_wrapper(simdata$Bhat,npc=npc,adjust = 'const',adj.const = s)



    list(data = simdata,
         mash.out=mash.out)

  },mc.cores = mc.cores)

  result

}



simu_mash_lfsr_udr = function(Ulist,pi,N,
                                 seed=12345,mean.range = 4,
                                 nreps = 20,mc.cores = 4,
                              adjust = NULL,s=0.1,
                                 prior='normal',df=10,half.uniform=FALSE){

  set.seed(seed)

  result = mclapply(1:nreps,function(i){
    # generate data
    simdata = simDataI.ult(N,Ulist,pi,prior=prior,
                           df=df,half.uniform=half.uniform,mean.range=mean.range)
    #datax = mash_set_data(simdata$Bhat,simdata$Shat)
    # fit model
    mash.out = mash_wrapper_udr(simdata$Bhat,adjust=adjust,adj.const =s)



    list(data = simdata,
         mash.out=mash.out)

  },mc.cores = mc.cores)

  result

}

get.fsr = function(out,alpha){

  fdr_result = lapply(out,function(x){
    #print(length(x))
    non_null_idx = which(x$data$B!=0)
    #print(non_null_idx)
    rej.mash = mashFDR(x$mash.out,alpha)
    fp.mash = c(fdp(rej.mash,non_null_idx),fsp(x$data$Bhat[rej.mash],x$data$B[rej.mash]),powr(rej.mash,non_null_idx))

    fp.mash
  })

  temp = do.call(rbind,fdr_result)

  colnames(temp) = c('FDR','FSR','power')

  temp
}


mash_wrapper_udr = function(Bhat,P=NULL,Shat=NULL,
                        nu = ncol(Bhat)+1,verbose=FALSE,
                        adjust = NULL,
                        U.true = NULL,
                        U.ed=NULL,
                        n.ed=NULL,
                        adj.const = 0.1){



  if(!is.null(P)){
    Z = sign(Bhat)*qnorm(1-P/2)
  }else{
    Z = Bhat
  }

  R = ncol(Z)
  N = nrow(Z)

  data = mash_set_data(Z,Shat)
  U.c = cov_canonical(data)

  if(is.null(U.true)){

    if(is.null(U.ed)){

      m.1by1 = mash_1by1(data)
      strong = get_significant_results(m.1by1)
      n.ed = length(strong)
      if(length(strong)==0){
        strong=1:N
        n.ed=N
      }
      #browser()
      init.random = list(U1 =drop(rWishart(1,R,diag(R))),
                         U2 =drop(rWishart(1,R,diag(R))),
                         U3 =drop(rWishart(1,R,diag(R))))
      U.init = ud_init(Bhat[strong,],n_rank1 = 1,U_unconstrained = init.random)
      U.ed = ud_fit(U.init,Bhat[strong,],verbose = FALSE)
    }
    U.est = U.ed$U
    if(!is.null(adjust)){

      if(adjust == 'prior'){
        U.est = lapply(1:length(U.est),function(k){
          nk = n.ed*U.ed$pi[k]
          if(nk>R){
            s2.hat = R*nu/(nk+nu)/sum(diag(solve(nk*(U.est[[k]]+diag(R)))))
            (nk/(nk+nu-R-1))*U.est[[k]] + ((R+1-nu+s2.hat)/(nk+nu-R-1))*diag(R)
          }else{
            NULL
          }

        })
      }

      if(adjust == "lb"){
        U.est = lapply(1:length(U.est),function(k){
          nk = n.ed*U.ed$pi[k]
          if(nk>1){
            U.est[[k]] + 2*diag(sqrt(2/nk),R)
            #U.est[[k]] + 2*diag(sqrt(2/nk*1/(diag(solve(U.est[[k]]+diag(R))))^2))
            #temp = eigen(U.est[[k]])
            #temp$vectors%*%(diag(pmax(2/sqrt(nk),temp$values)))%*%t(temp$vectors)
          }else{
            NULL
          }
        })
      }

      if(adjust == "const"){
        U.est = lapply(1:length(U.est),function(k){
          U.est[[k]] + diag(adj.const,R)
        })
      }

    }

    U.data = U.est[lengths(U.est) != 0]

  }else{
    U.data=U.true
  }

  out = mash(data,c(U.c,U.data),verbose=verbose)
  out

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

