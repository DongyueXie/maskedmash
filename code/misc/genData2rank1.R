

genData2 = function(nsamp = 100, err_sd=0.01,signal_sd = 2){
  ncond=5
  b1 = rnorm(nsamp,0,signal_sd)
  B.1 = matrix(cbind(b1,b1,0,0,0),nrow=nsamp, ncol=ncond) #independent effects
  b2 = rnorm(nsamp,0,signal_sd)
  B.2 = matrix(cbind(0,0,b2,b2,b2),nrow=nsamp, ncol=ncond) #independent effects
  
  B = rbind(B.1, B.2)
  
  Shat = matrix(err_sd, nrow=nrow(B), ncol=ncol(B))
  E = matrix(rnorm(length(Shat), mean=0, sd=Shat), nrow=nrow(B),ncol=ncol(B))
  Bhat = B+E
  row_ids = paste0("effect_", 1:nrow(B))
  col_ids = paste0("condition_", 1:ncol(B))
  rownames(B) = row_ids
  colnames(B) = col_ids
  rownames(Bhat) = row_ids
  colnames(Bhat) = col_ids
  rownames(Shat) = row_ids
  colnames(Shat) = col_ids
  return(list(B=B,Bhat=Bhat,Shat=Shat))
}