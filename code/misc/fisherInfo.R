n = 1000
U = matrix(c(2,1,1,2),nrow=2)
St = solve(U)
I = (kronecker(St,St) + tcrossprod(c(St)))/4
solve(n*I)

## mine calculation
St = solve(U)
Ii = matrix(0,nrow=3,ncol = 3)
Ii[1,1] = St[1,1]^2/2
Ii[1,2] = St[1,1]*St[1,2]
Ii[1,3] = St[1,2]^2/2
Ii[2,2] = St[1,2]^2 + St[1,1]*St[2,2]
Ii[2,3] = St[1,2]*St[2,2]
Ii[3,3] = St[2,2]^2/2
Ii = (Ii+t(Ii))-diag(diag(Ii))
solve(n*Ii)

set.seed(12345)
covs = array(dim = c(2,2,1000))
for(i in 1:1000){
  x = rmvnorm(n,sigma = U)
  covs[,,i] = cov(x)
}
var(covs[1,1,])
var(covs[2,2,])
var(covs[1,2,])
cov(covs[1,1,],covs[1,2,])
cov(covs[1,1,],covs[2,2,])
cov(covs[1,2,],covs[2,2,])

Ii2 = Ii
Ii2 = cbind(Ii2[,1:2],Ii2[,2],Ii2[,3])
Ii2 = rbind(Ii2[1:2,],Ii2[2,],Ii2[3,])
solve(Ii2)

Ii3 = matrix(0,nrow=4,ncol=4)
Ii3[1,1] = St[1,1]^2
Ii3[1,2] = St[1,1]*St[1,2]
Ii3[1,3] = St[1,1]*St[1,2]
Ii3[1,4] = St[1,2]^2
Ii3[2,2] = St[1,2]^2
Ii3[2,3] = St[1,1]*St[2,2]
Ii3[2,4] = St[1,2]*St[2,2]
Ii3[3,3] = St[1,2]^2
Ii3[3,4] = St[1,2]*St[2,2]
Ii3[4,4] = St[2,2]^2
Ii3 = (Ii3+t(Ii3))-diag(diag(Ii3))
Ii3 = Ii3/2
solve(n*Ii3)

