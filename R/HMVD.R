HMVD <-
function(Y,X,XX=NULL,weight = FALSE,C=1,nperm.max=0,model.type = 'C',method = 'test',adj=TRUE){
  ###############################check C
  if(!is.numeric(C)){stop('C can only be numeric !')}else{
    if(C<=0){stop('C should have a positive value!')}
  }
  #########check if Y is a vector
  if(!is.vector(Y)&&dim(Y)[2]!=1) stop('The outcome variable Y is not in vector form!')
  n = length(Y)
  if(model.type=='D'){if(range(Y)[1]<0|range(Y)[2]>1) stop('The outcome variable Y is not binary!')}
  ############check if the dimension of X is compatible with y
  if(is.null(X))stop('X is missing')
  if((is.vector(X)==FALSE) && (is.matrix(X)==FALSE))stop('X has to be a vector or a matrix')
  if(is.vector(X)){
    m=1;if(length(X)!=n) stop('X does not have the same length as the outcome variable Y!')
    ############ if there are missing covariates, use the mean for imputation
    id = which(is.na(X)); if(length(id)>0) X[id] = mean(X,na.rm=TRUE)
    X = as.matrix(X,ncol=1) ############### make X a matrix
  }else{
    if(is.matrix(X)){
      m = dim(X)[2];if(dim(X)[1]!=n) stop('The number of rows of X is not the length of the outcome variable Y!')
      for(i in 1:m){id = which(is.na(X[,i])); if(length(id)>0) X[id,i] = mean(X[,i],na.rm=TRUE)}
      v = c();for(i in 1:m) v[i]=var(X[,i])
      keep = which(v>0)
      if(length(keep)==0) stop('All columns of X are constant!')
      X_keep = matrix(X[,keep],nrow = n)
    }
  }
  pp = 0 ############## number of covariates
  if(!is.null(XX)){
    if((is.vector(XX)==FALSE) && (is.matrix(XX)==FALSE))stop('The covariates XX has to be a vector or a matrix!')
    ##############if XX is not NULL, then we need to check the dimension of XX as well
    if(is.vector(XX)){
      pp=1;if(length(XX)!=n) stop('The covariate XX does not have the same length as outcome variable!')
      id = which(is.na(XX)); if(length(id)>0) XX[id] = mean(XX,na.rm=TRUE)
    }else{
      if(is.matrix(XX)){
        pp = dim(XX)[2];if(dim(XX)[1]!=n) stop('The number of rows of XX is not the length of the outcome variable Y!')
        for(i in 1:pp){id = which(is.na(XX[,i])); if(length(id)>0) XX[id,i] = mean(XX[,i],na.rm=TRUE)}
      }
    }
  }
  
  #################check other paramters 
  if(!(model.type %in% c('C','D'))) stop('Wrong type of model type. Model.type can be either "C" (for continuous outcome) or "D" (for binary outcome).')
  #if(!(X.type %in% c('C','G'))) stop('Wrong type of X type. X.type can be either "C" (for continuous variable) or "G" (for SNPs).')
  if(length(weight)!=1|is.na(weight)) stop('Parameter "weight" can only be set to TRUE or FALSE.')
  if(!(method %in% c('test','estimation'))) stop('Wrong type of method. Model can be either "test" (for hypothesis testing) or "estimation" (for parameter estimation).')
  
  if(ncol(X_keep)>1){
    id2 = c()
    c = cor(X_keep);l=ncol(X_keep)
    for(i in 1:(l-1)){
      temp = c[(i+1):l,i]
      #### check whether variable (i+1):l is highly correlated with X_keep[,i]
      #### if yes, i.e. length(id)>0 remove all those variables but keep X_keep[,i]
      id = which(temp>.99)+i
      if(length(id)>0) id2 = c(id2,id[1])
    }
    id1 = setdiff(1:l,id2)
    keep2 = keep[id1]
    X_keep2 = matrix(X_keep[,id1],nrow=n)
  }
  
  if(length(keep2)==0) stop('All columns of X are removed because of high correlation!') ############## this should never happen
  
  X_final = X_keep2
  Y_final = Y
  ######### need the following steps if Y is binary.
  if(model.type=="D"){
    #################### first group snps if number of minor allele is less than 5. 
    #################### Maybe this step can be excluded in next version
    temp = X_keep2
    temp = matrix(temp,n)
    m = ncol(temp);mod = c();tmpcc = c()
    for(i in 1:m){tmp <- table(temp[,i]);tmpcc[i]=n-max(tmp)}

    while(sum(tmpcc<5)>0){
      po = which(tmpcc<5)[1]
      if(po<dim(temp)[2]) temp[,(po+1)] = temp[,po]+temp[,po+1]
      if(po==dim(temp)[2]) temp[,(po-1)] = temp[,po]+temp[,po-1]
      temp = temp[,-po]
      temp = matrix(temp,n)
      m = ncol(temp);mod=c();tmpcc=c()
      for(i in 1:m){tmp <- table(temp[,i]);tmpcc[i]=n-max(tmp)}
    }
    temp = matrix(temp,n)
    X_final = temp
    m = ncol(X_final)
    ################## we need to add two fake data points
    ################## such that glm never crashes.
    if(adj){
      Y_final = c(.5,.5,Y_final)
      X_final = rbind(apply(X_final,2,max),apply(X_final,2,min),X_final)
    }
  }

  out = list()

  #################### weight option is only applied to SNP variables
  ###################### update: weight option is now avaible for both types of variables.
  if(weight==TRUE){
    #if(SNP==FALSE) stop('Parameter "weight" can only be set to TRUE when X represents genotypes.')
    f = apply(X_final,2,sd)
    w = 1/sd
    X_final = t(t(X_final)*w)
  }

  m = ncol(X_final)
  if(m==0) stop('No variable is left in the group!')
  if(m==1){
    XXX = cbind(X_final,XX)
    if(model.type=='C'){fit = lm(Y_final~XXX); out$p.value = summary(fit)$coef[2,4]; out$theta = fit$coef[2]; out$theta_hat = fit$coef[2]}
    if(model.type=='D'){fit = glm(Y_final~XXX,family='quasibinomial'); out$p.value = summary(fit)$coef[2,4]; out$theta = fit$coef[2]; out$theta_hat = fit$coef[2]}
  }
  if(m>1){
    if(model.type=='C') out = HMM_lin(Y=Y_final,X = X_final,XX,C,nperm.max = nperm.max, method = method)
    if(model.type=='D') out = HMM_log(Y=Y_final,X = X_final,XX,C,nperm.max = nperm.max, method = method)  
  }
  out$keep1 = keep
  out$keep2 = keep2
  
  return(out)
}
