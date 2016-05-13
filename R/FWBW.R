FWBW <-
function(pi0_new,theta0_new,A0_new,x,mus,vars,m){
  a = matrix(0,m,2);
  b = matrix(0,m,2);b[m,]=1
  a[1,1] = pi0_new*dnorm(x[1],0,sqrt(vars[1])) #########pi0_new = P(Z1=0)
  a[1,2] = (1-pi0_new)*dnorm(x[1],theta0_new,sqrt(vars[1]))
  for(t in 2:m) for(i in 0:1) for(j in 0:1) 
    a[t,i+1] = a[t,i+1] + a[t-1,j+1]*A0_new[j+1,i+1]*dnorm(x[t],(i-j*mus[t])*theta0_new,sqrt(vars[t]))
  for(t in (m-1):1) for(i in 0:1) for(j in 0:1) 
   ## b[t,i+1] = b[t,i+1] + b[t+1,j+1]*A0_new[i+1,j+1]*dnorm(x[t+1],(i-j*mus[t+1])*theta0_new,sqrt(vars[t+1]))
    ## above is an error: should switch i and j
    b[t,i+1] = b[t,i+1] + b[t+1,j+1]*A0_new[i+1,j+1]*dnorm(x[t+1],(j-i*mus[t+1])*theta0_new,sqrt(vars[t+1]))
  return(cbind(a,b))
}
