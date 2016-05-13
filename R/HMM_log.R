HMM_log <-
function(Y,X,XX,C=1,nperm.max,method){
  n=length(Y);X = matrix(X,nrow=n)
  m=dim(X)[2];pp=0;if(!is.null(XX)) pp=dim(XX)[2]
 
  theta_hat = c(); out = list(); out$p.value=NA
  temp=matrix(0,n,m)
  #####theta_hat and its distribution
  ###S is the covariance matrix, and theta_hat are our data.
    theta_hat = c(); Res = matrix(0,n,m); Res1 = matrix(0,n,m); 
    V = list(); V1 = list(); V2 = list();
    sd_check = c()
    for(i in 1:m){
      if(is.null(XX)){XXX=X[,i]}else{XXX = cbind(X[,i],XX)}
      fit = glm(Y~XXX,family=quasibinomial)
      theta_hat[i] = fit$coef[2]
      sd_check[i] = summary(fit)$coef[2,2]
      temp[,i] = fit$fitted
      V[[i]] = cbind(rep(1,n),XXX)
      Res[,i] =  (Y-temp[,i])
      Res1[,i] = (1-temp[,i])
      V1[[i]] = V[[i]]*Res[,i]
      V2[[i]] = V[[i]]*temp[,i]*Res1[,i]
    }
    S = matrix(0,m,m);Sig = matrix(0,m,m)
    for(i in 1:m) for(j in 1:i){
      S1 = solve(t(V[[i]])%*%V2[[i]])
      S2 = solve(t(V[[j]])%*%V2[[j]])
      S3 = t(V1[[i]])%*%V1[[j]]
      S[i,j]=S[j,i] = (S1%*%(S3)%*%t(S2))[2,2]
    }
  
  ############check whether non-convergence happened....
  for(i in 1:m){
    flag = (abs(sd_check[i]-sqrt(S[i,i]))>.5*(abs(sd_check[i])))  
    if(flag) warning("glm.fit: algorithm did not converge. Suggest use adj=TRUE option")
  }
    
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
  power = -ceiling(log10(median(sqrt(vars))))-1
  power = 0 #### do not rescale theta_hat and vars
  ##### HMM_EM #####
  if(method=='test'){
    out = HMM_EM(theta_hat*10*power,mus,vars*100^power,x*10^power,const,C)
    if(nperm.max>0){
      nperm = 1/out$p.value*10; nperm = min(nperm,nperm.max)
      pperm = c()
      for(i in 1:nperm){
        theta_hat_perm = mvrnorm(1, mu=rep(0,m), S)
        x = c();x[1] = theta_hat_perm[1];x[2:m] = theta_hat_perm[2:m]-mus[2:m]*theta_hat_perm[2:m-1]
        pperm[i] = HMM_EM(theta_hat_perm*10^power,mus,vars*100^power,x*10^power,const,C)$p.value
      }
      #out$p.value.perm = mean(out$p.value>pperm)
      out$p.value.perm = mean(out$p.value>=c(pperm,out$p.value))
      out$nperm = nperm
    }
  }
  
  if(method=='estimation')  out = HMM_EM_fwbw(theta_hat,mus,vars,x,C) 
  
  return(out)
}
