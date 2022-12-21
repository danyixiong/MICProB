setwd("~/Dropbox/BioHPC/MICProB/Data_Code/Simulation")

#### functions ####
BMIR2_genX<-function(n,m,d,corr){
  #npoints<-n*m*d
  #x<-rnorm(npoints)
  #X<-matrix(x, ncol = d)
  
  # covariance matrix
  sigma = diag(1,nrow = d, ncol = d)
  sigma[upper.tri(sigma)] = corr
  sigma[lower.tri(sigma)] = corr
  
  X = mvtnorm::rmvnorm(n*m,sigma = sigma)
  probit_prob <- pnorm(b_true[1] + X %*% b_true[-1])
  prime_ind<-rbinom(n*m,1,probit_prob)
  #X_all<-list(X, prime_ind)
  
  bags<-list()
  j = 1
  for(i in 1:n){
    bags[[i]]<-list(X[j:(j+m-1),],prime_ind[j:(j+m-1)])
    j = j + m
  }
  return(bags)
}

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

#### check ####

# basic setting
# n, m, d, ppi
# n = 300; m=10; d=30; ppi=0.4

dd<-50 # number of replicates
mean.ppi<-rep(NA,dd)
mean.bag1<-rep(NA,dd)

set.seed(1779)

#b_true<-c(-3.5,rep(1,30)) # ppi=0.4
#beta_true<-c(0.5,rep(-1,15),rep(0.8,15))

b_true<-c(-1.4,rep(1,30))
beta_true<-c(0.5,rep(-1,15),rep(0.5,15))

b_true<-c(-1.5,rep(1,30)) # ppi=0.4
beta_true<-c(-2.5,rep(-1,15),rep(0.8,15))

b_true<-c(-3.5,rep(1,30)) # ppi=0.4
beta_true<-c(0.5,rep(-1,15),rep(0.8,15))

b_true<-c(-4,rep(1,30)) # ppi=0.4
beta_true<-c(0.5,rep(-1,15),rep(0.8,15))

# for rho = 0.95
b_true<-c(-2,rep(-1,15),rep(0.5,15)) # ppi=0.4
beta_true<-c(-1,rep(1,15),rep(-1,15))

for(i in 1:dd){
  test_X<-BMIR2_genX(n=300,m=10,d=length(b_true)-1,corr = 0.95)
  prime_prop<-unlist(lapply(test_X, function(x){mean(x[[2]])}))
  mean.ppi[i]<-mean(prime_prop)
  mean.bag1[i]<-mean(unlist(lapply(X=test_X, BMIR2_genY)))
}

summary(mean.ppi)
summary(mean.bag1)

#### generate data ####
nrep<-50 # number of replicates
save_f<-"revision/sim_corr_2/"
corr_vec = c(0,0.2,0.4,0.6)

#### training bags ####

set.seed(1779)

for(i in 1:length(corr_vec)){
  
  mean.ppi<-rep(NA,nrep)
  mean.bag1<-rep(NA,nrep)
  
  if(i==1){
    b_true<-c(-1.5,rep(1,30)) # ppi=0.4
    beta_true<-c(-2.5,rep(-1,15),rep(0.8,15))
  } else if(i==2){
    b_true<-c(-3.5,rep(1,30)) # ppi=0.4
    beta_true<-c(0.5,rep(-1,15),rep(0.8,15))
  } else if(i==3){
    b_true<-c(-4,rep(1,30)) # ppi=0.4
    beta_true<-c(0.5,rep(-1,15),rep(0.8,15))
  } else{
    b_true<-c(-3.5,rep(1,30)) # ppi=0.4
    beta_true<-c(0.5,rep(-1,15),rep(0.8,15))
  }
  
  for(j in 1:nrep){
    #X_all<-BMIR2_genX(n=300,m=10,d=30,corr=corr_vec[i])
    X_all<-BMIR2_genX(n=300,m=10,d=30,corr=0.01)
    X<-lapply(X_all, function(x){x[[1]]})
    delta_true<-lapply(X_all, function(x){x[[2]]})
    
    y = unlist(lapply(X_all, BMIR2_genY))
    
    prime_prop<-unlist(lapply(X_all, function(x){mean(x[[2]])}))
    mean.ppi[j]<-mean(prime_prop)
    mean.bag1[j]<-mean(y)
    
    #save(X,y,beta_true,b_true, delta_true,file = paste0(save_f,"/","train_",5,"_",j, ".Rdata"))
  }
  
  cat("Scenario:",i,"created. ",mean(mean.ppi),';',mean(mean.bag1),"\n")
}

#### test bags ####

for(i in 1:length(corr_vec)){
  
  mean.ppi<-rep(NA,nrep)
  mean.bag1<-rep(NA,nrep)
  
  for(j in 1:nrep){
    #X_all<-BMIR2_genX(n=300,m=10,d=30,corr=corr_vec[i])
    X_all<-BMIR2_genX(n=300,m=10,d=30,corr=0.01)
    X<-lapply(X_all, function(x){x[[1]]})
    delta_true<-lapply(X_all, function(x){x[[2]]})
    
    y = unlist(lapply(X_all, BMIR2_genY))
    
    prime_prop<-unlist(lapply(X_all, function(x){mean(x[[2]])}))
    mean.ppi[j]<-mean(prime_prop)
    mean.bag1[j]<-mean(y)
    
    #save(X,y,beta_true,b_true, delta_true,file = paste0(save_f,"/","test_",5,"_",j, ".Rdata"))
  }
  
  cat("Scenario:",i,"created. ",mean(mean.ppi),';',mean(mean.bag1),"\n")
}

load("revision/sim_corr/test_5_1.Rdata")
mean(do.call(rbind, delta_true))
mean(y)

load("revision/sim_corr/train_5_1.Rdata")
mean(do.call(rbind, delta_true))
mean(y)

load("revision/sim_corr/train_5_1.Rdata")
X = do.call(rbind, X)
cor(X)

#### rho = 0.95 ####
save_f = "revision/sim_corr_2/"

b_true<-c(-2,rep(-1,15),rep(0.5,15)) # ppi=0.39
beta_true<-c(-1,rep(1,15),rep(-1,15)) # p = 0.27

set.seed(1779)

mean.ppi<-rep(NA,nrep)
mean.bag1<-rep(NA,nrep)

for(j in 1:nrep){
  #X_all<-BMIR2_genX(n=300,m=10,d=30,corr=corr_vec[i])
  X_all<-BMIR2_genX(n=300,m=10,d=30,corr=0.95)
  X<-lapply(X_all, function(x){x[[1]]})
  delta_true<-lapply(X_all, function(x){x[[2]]})
  
  y = unlist(lapply(X_all, BMIR2_genY))
  
  prime_prop<-unlist(lapply(X_all, function(x){mean(x[[2]])}))
  mean.ppi[j]<-mean(prime_prop)
  mean.bag1[j]<-mean(y)
  
  save(X,y,beta_true,b_true, delta_true,file = paste0(save_f,"/","test_",6,"_",j, ".Rdata"))
}

mean(mean.ppi)
mean(mean.bag1)

load("revision/sim_corr_2/train_6_1.Rdata")
load("revision/sim_corr_2/test_6_1.Rdata")


