

#'@param X data matrix N*R
#'@param w posterior weight matrix N*K
#'@param U a list of estimated U
#'@param V an array of dimension (R,R,N)
#'@return info mat and its inverse of each estimated U
calc_info_mat = function(X,w,U,V){

  K = length(U)
  N = nrow(X)
  R = ncol(X)

  if(length(dim(V))==2){
    V = array(V,dim=c(R,R,N))
  }

  U.var = list()

  for(k in 1:K){

    S_tilde_k = V
    for(i in 1:N){
      S_tilde_k[,,i] = solve(V[,,i]+U[[k]])
    }

    A = c()


    for(r in 1:R){
      for(jr in r:(R)){

        a = c()


        for(l in 1:R){
          for(kl in l:(R)){

            if(r==jr){
              if(l==kl){

                a = c(a, sum(w[,k]*(S_tilde_k[r,l,])^2)/2)


              }else{

                a = c(a, sum(w[,k]*S_tilde_k[r,l,]*S_tilde_k[r,kl,]))


              }
            }else{
              if(l==kl){

                a = c(a, sum(w[,k]*S_tilde_k[l,r,]*S_tilde_k[l,jr,]))


              }else{
                a = c(a, sum(w[,k]*(S_tilde_k[l,jr,]*S_tilde_k[r,kl,] + S_tilde_k[r,l,]*S_tilde_k[jr,kl,])))

              }
            }

          }
        }

        A = rbind(A,a)

      }
    }





    temp = matrix(0,nrow=R,ncol=R)
    temp[lower.tri(temp,diag = T)] = diag(solve(A))
    temp = (temp+t(temp))-diag(diag(temp))
    U.var[[k]] = temp


  }

  U.var

}
