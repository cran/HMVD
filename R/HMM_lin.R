HMM_lin <-
function(Y,X,XX,weight = TRUE,C=1,nperm.max,method){
  n=length(Y);m=dim(X)[2];pp=0
  if(!is.null(XX)) pp=dim(XX)[2]
############## shift X and Y to mean 0 only for linear regression.
  X = scale(X,scale=FALSE)
  Y = Y-mean(Y)

  theta_hat = c(); 
  Res = matrix(0,n,m); 
  for(i in 1:m){
    if(is.null(XX)){XXX = X[,i]}else{XXX = cbind(X[,i],XX)}
    th = solve(t(XXX)%*%XXX)%*%(t(XXX)%*%Y)
    theta_hat[i] = th[1]
    Res[,i] = Y-XXX%*%th
  }
  
  S = matrix(0,m,m);Sig = matrix(0,m,m)
  for(i in 1:m) for(j in 1:i){
    if(is.null(XX)){
      XXX1=X[,i];XXX2=X[,j]
      S[i,j]=S[j,i]=(1/sum(XXX1^2)*sum(XXX1*XXX2)*1/(sum(XXX2^2)))
    }else{
      XXX1=cbind(X[,i],XX);XXX2=cbind(X[,j],XX)
      S[i,j]=S[j,i]=(solve(t(XXX1)%*%XXX1)%*%(t(XXX1)%*%XXX2)%*%solve(t(XXX2)%*%XXX2))[1,1]
    }
    Sig[i,j] = Sig[j,i] = sum(Res[,i]*Res[,j])/n
  }
  S=S*Sig
  
  ###conditional distribution
  ###step 1: find the conditional mean and variance mus = s[t,p]/s[p,p]
  mus = c(); mus[1] = 0; for(i in 2:m) mus[i] = S[i,i-1]/S[i-1,i-1]
  vars = c(); vars[1] = S[1,1]; for(i in 2:m) vars[i] = S[i,i] - S[i,i-1]^2/S[i-1,i-1]
  x = c();x[1] = theta_hat[1]
  x[2:m] = theta_hat[2:m]-mus[2:m]*theta_hat[2:m-1]
  
  SS1 = S
  SS2 = t(t(cbind(rep(0,m),S[,-m]))*mus)
  SS3 = S; SS3[1,]=SS3[,1] = 0; SS3[2:m,2:m] = S[-m,-m]; SS3 = SS3*(mus%*%t(mus))
  Omega = SS1-SS2-t(SS2)+SS3
  Sig = matrix(0,m,m)
  for(i in 1:m) Sig[i,] = Omega[i,]/vars[i]/vars
  Sig = Sig*((1-mus)%*%t(1-mus))
  const = sum(diag(Sig))/sum(Sig)
  
  ############# determine how to rescale the thata_hat and the vars
  #power = -ceiling(log10(median(sqrt(vars))))-1
  power = 0 #### do not rescale theta_hat and vars

  ##### HMM_EM #####
  if(method=='test'){
    out = HMM_EM(theta_hat*10^power,mus,vars*100^power,x*10^power,const,C)
    if(nperm.max>0){
      nperm = 1/out$p.value*10; nperm = min(nperm,nperm.max)
      pperm = c()
      for(i in 1:nperm){
        theta_hat_perm = mvrnorm(1, mu=rep(0,m), S)
        x = c();x[1] = theta_hat_perm[1];x[2:m] = theta_hat_perm[2:m]-mus[2:m]*theta_hat_perm[2:m-1]
        pperm[i] = HMM_EM(theta_hat_perm*10^power,mus,vars*100^power,x*10^power,const,C)$p.value
      }
#      out$p.value.perm = mean(out$p.value>pperm)
      out$p.value.perm = mean(out$p.value>=c(pperm,out$p.value))
      out$nperm = nperm
    }
  }

  if(method=='estimation')  out = HMM_EM_fwbw(theta_hat,mus,vars,x,C) 
  
  return(out)
}
