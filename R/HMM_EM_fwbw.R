HMM_EM_fwbw <-
function(theta_hat,mus,vars,x,C){
  ##### HMM #####
  ### initialize: 
  m = length(theta_hat)
  pi0 = .5;A0 = matrix(c(0.5,.5,.5,.5),2,2)
  t1 = mean(sort(theta_hat)[1:max(m*.2,1)])
  t2 = mean(sort(theta_hat)[m:(m*.8)])
  if((t1+t2)>0){theta_st = t2}else{theta_st = t1}
  pi0_new = pi0; A0_new=A0; theta0_new=theta_st; A0_old=matrix(0,2,2);theta0_old=-10
  iter=0

  ### M-step:
  while((abs(theta0_new-theta0_old))>1e-10&&(iter<1000)){
    iter=iter+1;
    pi0_old = pi0_new; A0_old=A0_new; theta0_old=theta0_new
    #print(iter)
    #############update theta
    R = FWBW(pi0_new,theta0_new,A0_new,x,mus,vars,m); a = R[,1:2];b=R[,3:4]
    Pm2 = matrix(0,m-1,4)
    Pm1 = matrix(0,m,2)
    for(t in 2:m) for(i in 0:1) for(j in 0:1)
      Pm2[t-1,2*i+j+1] = a[t-1,i+1]*A0_new[i+1,j+1]*
      dnorm(x[t],(j-mus[t]*i)*theta0_new,sqrt(vars[t]))*b[t,j+1]
    for(i in 0:1)
      Pm1[1,i+1] = a[1,i+1]*b[1,i+1] 
    for(t in 2:m){Pm1[t,1]=Pm2[t-1,1]+Pm2[t-1,3];Pm1[t,2]=Pm2[t-1,2]+Pm2[t-1,4]}
#log(sum(Pm1[1,]))+sum(log(A0_new[,2]))
   
    theta_upp = sum((Pm2[,2]-Pm2[,3]*mus[2:m]+Pm2[,4]*(1-mus[2:m]))*(theta_hat[2:m]-mus[2:m]*theta_hat[1:(m-1)])/vars[2:m]);
    theta_upp = theta_upp + Pm1[1,2]*theta_hat[1]/vars[1]
    theta_low = sum((Pm2[,2]+Pm2[,3]*mus[2:m]^2+Pm2[,4]*(1-mus[2:m])^2)/vars[2:m]);
    theta_low = theta_low + Pm1[1,2]/vars[1]
    theta0_new = theta_upp/theta_low
    
###################update pi
    R = FWBW(pi0_new,theta0_new,A0_new,theta_hat,mus,vars,m); a = R[,1:2];b=R[,3:4]
    for(t in 2:m) for(i in 0:1) for(j in 0:1)
      Pm2[t-1,2*i+j+1] = a[t-1,i+1]*A0_new[i+1,j+1]*
      dnorm(x[t],(j-mus[t]*i)*theta0_new,sqrt(vars[t]))*b[t,j+1]
    for(i in 0:1)
      Pm1[1,i+1] = a[1,i+1]*b[1,i+1] 
    for(t in 2:m){Pm1[t,1]=Pm2[t-1,1]+Pm2[t-1,3];Pm1[t,2]=Pm2[t-1,2]+Pm2[t-1,4]}
    pi0_new = (Pm1[1,1])/(Pm1[1,2]+Pm1[1,1])


##################update A
    R = FWBW(pi0_new,theta0_new,A0_new,theta_hat,mus,vars,m); a = R[,1:2];b=R[,3:4]
    for(t in 2:m) for(i in 0:1) for(j in 0:1)
      Pm2[t-1,2*i+j+1] = a[t-1,i+1]*A0_new[i+1,j+1]*dnorm(x[t],(j-mus[t]*i)*theta0_new,sqrt(vars[t]))*b[t,j+1]
    for(i in 0:1)
        Pm1[1,i+1] = a[1,i+1]*b[1,i+1] 
    for(t in 2:m){Pm1[t,1]=Pm2[t-1,1]+Pm2[t-1,3];Pm1[t,2]=Pm2[t-1,2]+Pm2[t-1,4]}
    #log(sum(Pm1[1,]))+sum(log(A0_new[,2]))

    A0_new[1,1] = sum(Pm2[,1])/(sum(Pm2[,1])+sum(Pm2[,2])+sum(Pm2[1,]))
    A0_new[2,1] = sum(Pm2[,3])/(sum(Pm2[,3])+sum(Pm2[,4])+sum(Pm2[1,]))    
    A0_new[,2] = 1-A0_new[,1]

    if(is.na(theta0_new)){break}
  }
  PPost=Pm1[,2]/(Pm1[,1]+Pm1[,2])
  Ptheta=theta0_new
FF = matrix(0,4,m)
FF[1,] = dnorm(x,0,sqrt(vars))
FF[2,] = dnorm(x,theta0_old,sqrt(vars))
FF[3,] = dnorm(x,-mus*theta0_old,sqrt(vars))
FF[4,] = dnorm(x,(1-mus)*theta0_old,sqrt(vars))
Ppi = A0_old[2,1]/(A0_old[2,1]+A0_old[1,2])
a = c(Ppi*A0_old[1,1],Ppi*A0_old[1,2],(1-Ppi)*A0_old[2,1],(1-Ppi)*A0_old[2,2])
ln1 = sum(log(colSums(FF*a)))+sum(log(A0_old[,2]))*C

  out1 = list()
  out1$prob = PPost
  out1$A = A0_new
  out1$theta = Ptheta
  out1$iter = iter

  
  ##### HMM #####
  ### initialize: 
  pi0 = .5;A0 = matrix(c(0.0001,0.0001,0.9999,0.9999),2,2);theta0 = 0;
  t1 = mean(sort(theta_hat)[1:max(m*.2,1)])
  t2 = mean(sort(theta_hat)[m:(m*.8)])
  if((t1+t2)>0){theta_st = t2}else{theta_st = t1}
  pi0_new = pi0; A0_new=A0; theta0_new=theta_st; A0_old=matrix(0,2,2);theta0_old=-10
  iter=0
  ### M-step:
  while((abs(theta0_new-theta0_old))>1e-10&&(iter<1000)){
    iter=iter+1;
    ########## update theta0
    A0_old = A0_new;theta0_old = theta0_new;pi0_old = pi0_new;
    #print(iter)
    R = FWBW(pi0_new,theta0_new,A0_new,theta_hat,mus,vars,m); a = R[,1:2];b=R[,3:4]
    Pm2 = matrix(0,m-1,4)
    Pm1 = matrix(0,m,2)
    for(t in 2:m) for(i in 0:1) for(j in 0:1)
      Pm2[t-1,2*i+j+1] = a[t-1,i+1]*A0_new[i+1,j+1]*
      dnorm(x[t],(j-mus[t]*i)*theta0_new,sqrt(vars[t]))*b[t,j+1]
    for(i in 0:1)
      Pm1[1,i+1] = a[1,i+1]*b[1,i+1] 
    for(t in 2:m){Pm1[t,1]=Pm2[t-1,1]+Pm2[t-1,3];Pm1[t,2]=Pm2[t-1,2]+Pm2[t-1,4]}
    theta_upp = sum((Pm2[,2]-Pm2[,3]*mus[2:m]+Pm2[,4]*(1-mus[2:m]))*(theta_hat[2:m]-mus[2:m]*theta_hat[1:(m-1)])/vars[2:m]);
    theta_upp = theta_upp + Pm1[1,2]*theta_hat[1]/vars[1]
    theta_low = sum((Pm2[,2]+Pm2[,3]*mus[2:m]^2+Pm2[,4]*(1-mus[2:m])^2)/vars[2:m]);
    theta_low = theta_low + Pm1[1,2]/vars[1]
    theta0_new = theta_upp/theta_low
########## update pi0 and A0
    R = FWBW(pi0_new,theta0_new,A0_new,theta_hat,mus,vars,m); a = R[,1:2];b=R[,3:4]
    Pm2 = matrix(0,m-1,4)
    Pm1 = matrix(0,m,2)
    for(t in 2:m) for(i in 0:1) for(j in 0:1)
      Pm2[t-1,2*i+j+1] = a[t-1,i+1]*A0_new[i+1,j+1]*
      dnorm(x[t],(j-mus[t]*i)*theta0_new,sqrt(vars[t]))*b[t,j+1]
    for(i in 0:1) Pm1[1,i+1] = a[1,i+1]*b[1,i+1] 
    for(t in 2:m){Pm1[t,1]=Pm2[t-1,1]+Pm2[t-1,3];Pm1[t,2]=Pm2[t-1,2]+Pm2[t-1,4]}

    pi0_new = (Pm1[1,1])/(Pm1[1,2]+Pm1[1,1])
    A0_new[1,1] = sum(Pm2[,1])/(sum(Pm2[,1])+sum(Pm2[,2])+sum(Pm2[1,]))
    A0_new[2,1] = sum(Pm2[,3])/(sum(Pm2[,3])+sum(Pm2[,4])+sum(Pm2[1,]))    
    A0_new[,2] = 1-A0_new[,1]

    if(is.na(theta0_new)){break}
  }
  PPost=Pm1[,2]/(Pm1[,1]+Pm1[,2])
  Ptheta=theta0_new
FF = matrix(0,4,m)
FF[1,] = dnorm(x,0,sqrt(vars))
FF[2,] = dnorm(x,theta0_old,sqrt(vars))
FF[3,] = dnorm(x,-mus*theta0_old,sqrt(vars))
FF[4,] = dnorm(x,(1-mus)*theta0_old,sqrt(vars))
Ppi = A0_old[2,1]/(A0_old[2,1]+A0_old[1,2])
a = c(Ppi*A0_old[1,1],Ppi*A0_old[1,2],(1-Ppi)*A0_old[2,1],(1-Ppi)*A0_old[2,2])
ln2 = sum(log(colSums(FF*a)))+sum(log(A0_old[,2]))*C

  out2 = list()
  out2$prob = PPost
  out2$theta = Ptheta
  out2$A = A0_new
  out2$iter = iter


  if(is.na(ln1)) return(out2)
  if(is.na(ln2)) return(out1)
  if(ln1>ln2){return(out1)}else{return(out2)}


  #out = list()
  #out$out1 = out1; out$out2 = out2; out$ln1 = ln1; out$ln2 = ln2
  #return(out)
}
