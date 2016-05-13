HMM_EM <-
function(theta_hat,mus,vars,x,const,C){
  ##### HMM_EM #####
  m = length(theta_hat)
  t1 = mean(sort(theta_hat)[1:max(m*.2,1)])
  t2 = mean(sort(theta_hat)[m:(m*.8)])
  if((t1+t2)>0){theta0_st = t2}else{theta0_st = t1}
  
  W00_log = dnorm(x,0,sqrt(vars),log=TRUE)
  l_d = sum(W00_log)
#  W00 = dnorm(x,0,sqrt(vars))
  #############FOR OUT1
  ### initialize: 
  theta0_old = -999;a01_old = .99;a10_old = .99
  theta0_new = theta0_st;a01_new = .5;a10_new = .5
  ### M step
  ct1 = 0
  while(sum(abs(theta0_new-theta0_old))>1e-5 & ct1<1000){
    ct1 = ct1+1
    theta0_old = theta0_new; a01_old = a01_new; a10_old = a10_new;
    a00_old = 1-a01_old; a11_old = 1-a10_old
    p_old = a01_old/(a01_old+a10_old)
    
    pi00 = (1-p_old)*(1-a01_old)  
    pi01 = (1-p_old)*a01_old
    pi10 = p_old*a10_old
    pi11 = p_old*(1-a10_old)
    W01_log = dnorm(x,theta0_old,sqrt(vars),log=TRUE)
    W10_log = dnorm(x,-mus*theta0_old,sqrt(vars),log=TRUE)
    W11_log = dnorm(x,(1-mus)*theta0_old,sqrt(vars),log=TRUE) ##works for the first one since mus[1]=0
    W_adj = cbind(W00_log,W01_log,W10_log,W11_log); WA = apply(W_adj,1,max)
    W_adj = W_adj-WA
    W = cbind(exp(W_adj[,1])*pi00,exp(W_adj[,2])*pi01,exp(W_adj[,3])*pi10,exp(W_adj[,4])*pi11);W = W/apply(W,1,sum)
    
    #### updata theta0_new
    upp = (W[,2]-W[,3]*mus+W[,4]*(1-mus))*x/vars
    low = (W[,2]+W[,3]*mus^2+W[,4]*(1-mus)^2)/vars
    theta0_new = sum(upp)/sum(low)
    
    #### updata transition matrix a10 and a01
    W01_log = dnorm(x,theta0_new,sqrt(vars),log=TRUE)
    W10_log = dnorm(x,-mus*theta0_new,sqrt(vars),log=TRUE)
    W11_log = dnorm(x,(1-mus)*theta0_new,sqrt(vars),log=TRUE) ##works for the first one since mus[1]=0
    W_adj = cbind(W00_log,W01_log,W10_log,W11_log); WA = apply(W_adj,1,max)
    W_adj = W_adj-WA
    W = cbind(exp(W_adj[,1])*pi00,exp(W_adj[,2])*pi01,exp(W_adj[,3])*pi10,exp(W_adj[,4])*pi11);W = W/apply(W,1,sum)
    A = sum(W[,1]); B = sum(W[,2])+sum(W[,3]); D = sum(W[,4])
    y = a10_old;z = a01_old
    for(iii in 1:3){
      ##### update z
      L1 = (A+B+C+D)*y+(A-C)
      z_new = (-L1+sqrt(L1^2+4*C*(B+C+D)*y))/2/C
      #C*z_new^2 + L1*z_new-(B+C+D)*y;z_new - z
      z = z_new
      #################### update y
      L2 = (A+B+C+D)*z+D
      y_new = (-L2+sqrt(L2^2+4*C*(A+B)*z))/2/C
      #y_new - y
      y = y_new
    }
    a01_new = z
    a10_new = y
    a00_new=1-a01_new; a11_new = 1-a10_new
  }
  W01_log = dnorm(x,theta0_old,sqrt(vars),log=TRUE)
  W10_log = dnorm(x,-mus*theta0_old,sqrt(vars),log=TRUE)
  W11_log = dnorm(x,(1-mus)*theta0_old,sqrt(vars),log=TRUE) ##works for the first one since mus[1]=0
  W_adj = cbind(W00_log,W01_log,W10_log,W11_log); WA = apply(W_adj,1,max)
  W_adj = W_adj-WA
  W = cbind(exp(W_adj[,1])*pi00,exp(W_adj[,2])*pi01,exp(W_adj[,3])*pi10,exp(W_adj[,4])*pi11)
  P = apply(W,1,sum)
  W = W/P

  p_new = a01_new/(a01_new+a10_new)
  pi00 = (1-p_new)*(1-a01_new)  
  pi01 = (1-p_new)*a01_new
  pi10 = p_new*a10_new
  pi11 = p_new*(1-a10_new)

  l_n = sum(log(P))+sum(WA)+C*log(a01_new*a11_new)
  LRT = (2*l_n-2*l_d)*const
  Pvalue = pchisq(LRT,df=1,lower.tail=FALSE)
  out1 = list()
  out1$p.value = Pvalue
  
  ####################### now start from a different start point #################
  #######################  FOR OUT2 ###################
  theta0_old = -999;a01_old = .5;a10_old = .5
  theta0_new = theta0_st;a01_new = .9999;a10_new = .0001
  ### M step
  ct2 = 0
  while(sum(abs(theta0_new-theta0_old))>1e-7 & ct2<1000){
    ct2 = ct2+1
    theta0_old = theta0_new; a01_old = a01_new; a10_old = a10_new;
    a00_old = 1-a01_old; a11_old = 1-a10_old
    p_old = a01_old/(a01_old+a10_old)
    
    pi00 = (1-p_old)*(1-a01_old)  
    pi01 = (1-p_old)*a01_old
    pi10 = p_old*a10_old
    pi11 = p_old*(1-a10_old)
    W01_log = dnorm(x,theta0_old,sqrt(vars),log=TRUE)
    W10_log = dnorm(x,-mus*theta0_old,sqrt(vars),log=TRUE)
    W11_log = dnorm(x,(1-mus)*theta0_old,sqrt(vars),log=TRUE) ##works for the first one since mus[1]=0
    W_adj = cbind(W00_log,W01_log,W10_log,W11_log); WA = apply(W_adj,1,max)
    W_adj = W_adj-WA
    W = cbind(exp(W_adj[,1])*pi00,exp(W_adj[,2])*pi01,exp(W_adj[,3])*pi10,exp(W_adj[,4])*pi11)
    P = apply(W,1,sum)
    W = W/P

    ################### update theta0_new #############################
    upp = (W[,2]-W[,3]*mus+W[,4]*(1-mus))*x/vars
    low = (W[,2]+W[,3]*mus^2+W[,4]*(1-mus)^2)/vars
    theta0_new = sum(upp)/sum(low)
    
    ################### update transition matrix ##########################
    W01_log = dnorm(x,theta0_new,sqrt(vars),log=TRUE)
    W10_log = dnorm(x,-mus*theta0_new,sqrt(vars),log=TRUE)
    W11_log = dnorm(x,(1-mus)*theta0_new,sqrt(vars),log=TRUE) ##works for the first one since mus[1]=0
    W_adj = cbind(W00_log,W01_log,W10_log,W11_log); WA = apply(W_adj,1,max)
    W_adj = W_adj-WA
    W = cbind(exp(W_adj[,1])*pi00,exp(W_adj[,2])*pi01,exp(W_adj[,3])*pi10,exp(W_adj[,4])*pi11)
    P = apply(W,1,sum)
    W = W/P
    A = sum(W[,1]); B = sum(W[,2])+sum(W[,3]); D = sum(W[,4])
    y = a10_old;z = a01_old
    
    for(iii in 1:3){
      ##### update z
      L1 = (A+B+C+D)*y+(A-C)
      z_new = (-L1+sqrt(L1^2+4*C*(B+C+D)*y))/2/C
      #C*z_new^2 + L1*z_new-(B+C+D)*y;z_new - z
      z = z_new
      #################### update y
      L2 = (A+B+C+D)*z+D
      y_new = (-L2+sqrt(L2^2+4*C*(A+B)*z))/2/C;
      #y_new - y
      y = y_new
    }
    a01_new = z
    a10_new = y
    a00_new=1-a01_new; a11_new = 1-a10_new
  }
  W01_log = dnorm(x,theta0_old,sqrt(vars),log=TRUE)
  W10_log = dnorm(x,-mus*theta0_old,sqrt(vars),log=TRUE)
  W11_log = dnorm(x,(1-mus)*theta0_old,sqrt(vars),log=TRUE) ##works for the first one since mus[1]=0
  W_adj = cbind(W00_log,W01_log,W10_log,W11_log); WA = apply(W_adj,1,max)
  W_adj = W_adj-WA
  W = cbind(exp(W_adj[,1])*pi00,exp(W_adj[,2])*pi01,exp(W_adj[,3])*pi10,exp(W_adj[,4])*pi11)
  P = apply(W,1,sum)
  W = W/P

  p_new = a01_new/(a01_new+a10_new)
  pi00 = (1-p_new)*(1-a01_new)  
  pi01 = (1-p_new)*a01_new
  pi10 = p_new*a10_new
  pi11 = p_new*(1-a10_new)
  l_n = sum(log(P))+sum(WA)+C*log(a01_new*a11_new) #########WA is added b/c P is nolonger the true probability 
  LRT = (2*l_n-2*l_d)*const
  Pvalue = pchisq(LRT,df=1,lower.tail=FALSE)
  out2 = list()
  out2$p.value = Pvalue

  #### output the solution that has smaller pvalue/ better convergence
  if(out2$p.value<out1$p.value){return(out2)}else{return(out1)}
}
