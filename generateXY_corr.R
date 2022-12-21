## generate multiple instance data with correlated features

#### generate instances ####
BMIR2_genX<-function(n,m,d,corr){
  # covariance matrix
  sigma = diag(1,nrow = d, ncol = d)
  sigma[upper.tri(sigma)] = corr
  sigma[lower.tri(sigma)] = corr
  
  X = mvtnorm::rmvnorm(n*m,sigma = sigma)
  probit_prob <- pnorm(b_true[1] + X %*% b_true[-1])
  prime_ind<-rbinom(n*m,1,probit_prob)
  
  bags<-list()
  j = 1
  for(i in 1:n){
    bags[[i]]<-list(X[j:(j+m-1),],prime_ind[j:(j+m-1)])
    j = j + m
  }
  return(bags)
}

#### generate bag label ####
BMIR2_genY<-function(X_all){
  # covariate matrix
  X<-X_all[[1]]
  
  # primary indicator
  prime_ind<-X_all[[2]]
  
  if(sum(prime_ind)==0){ # current bag has no primary instance
    X.prime<-rep(0,ncol(X))
  } else{ # current bag has at least one primary instance
    # design matrix formed by primary instances
    X.prime<-X[prime_ind==1,,drop=FALSE]
  }
  
  z_mu<-beta_true[1] + sum(X.prime %*% beta_true[-1])
  
  # probability of positive bag
  probit_prob<-pnorm(z_mu)
  y<-rbinom(1,1,probit_prob)
  
  return(y)
}

#### simulated data ####
nrep = 50 # number of replicates
corr_vec = c(0.2,0.4,0.6,0.95)

for(i in 1:length(corr_vec)){
  
  mean.ppi<-rep(NA,nrep)
  mean.bag1<-rep(NA,nrep)
  
  if(i %in% 1:3){
    b_true<-c(-3.5,rep(1,30))
    beta_true<-c(0.5,rep(-1,15),rep(0.8,15))
  } else{
    b_true<-c(-2,rep(-1,15),rep(0.5,15))
    beta_true<-c(-1,rep(1,15),rep(-1,15))
  }
  
  for(j in 1:nrep){
    X_all<-BMIR2_genX(n=300,m=10,d=30,corr=corr_vec[i])
    X<-lapply(X_all, function(x){x[[1]]})
    delta_true<-lapply(X_all, function(x){x[[2]]})
    
    y = unlist(lapply(X_all, BMIR2_genY))
    
    prime_prop<-unlist(lapply(X_all, function(x){mean(x[[2]])}))
    mean.ppi[j]<-mean(prime_prop)
    mean.bag1[j]<-mean(y)
    
    save(X,y,beta_true,b_true, delta_true,file = paste0("simdat_corr_",i,"_",j, ".Rdata"))
  }
  
  cat("Scenario:",i,"created. ",mean(mean.ppi),';',mean(mean.bag1),"\n")
}

