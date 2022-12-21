## generate multiple instance data with independent features

#### generate instances ####
BMIR2_genX<-function(n,m,d){
  npoints<-n*m*d
  x<-rnorm(npoints)
  X<-matrix(x, ncol = d)
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
sim_param<-matrix(c(
  # n = 150, 300, 450, 600
  150,10,30,0.4, # 1
  300,10,30,0.4, #2 basic setting
  450,10,30,0.4, #3
  600,10,30,0.4, #4
  # m = 5, 10, 20, 40
  300,5,30,0.4, #5
  300,10,30,0.4, #6 basic setting
  300,20,30,0.4, #7
  300,40,30,0.4, #8
  # d = 2, 15, 30, 45
  300,10,2,0.4, #9
  300,10,15,0.4,#10
  300,10,30,0.4, #11 basic setting
  300,10,45,0.4, #12
  # ppi = 0.1, 0.4, 0.6, 0.9
  300,10,30,0.1, #13
  300,10,30,0.4, #14 basic setting
  300,10,30,0.6, #15
  300,10,30,0.9) #16
  ,byrow=T,ncol=4,
  dimnames = list(NULL,c("n","m","d","ppi")))

nrep<-50 # number of replicates

for(i in 1:nrow(sim_param)){

  if(i %in% c(1:8,11,14)){
    b_true<-c(-1.4,rep(1,30))
    beta_true<-c(0.5,rep(-1,15),rep(0.5,15))
    
  } else if(i==9){
    b_true<-c(-0.4,rep(1,2))
    beta_true<-c(0.5,rep(-1,1),rep(0.5,1)) 
    
  } else if(i==10){
    b_true<-c(-1,rep(1,15)) 
    beta_true<-c(0.5,rep(-1,8),rep(0.5,7)) 
    
  } else if(i==12){
    b_true<-c(-1.7,rep(1,45)) 
    beta_true<-c(0.5,rep(-1,23),rep(0.5,22)) 
    
  } else if(i==13){
    b_true<-c(-7,rep(1,30)) 
    beta_true<-c(0.5,rep(-1,15),rep(0.5,15)) 
    
  } else if(i==15){
    b_true<-c(1.4,rep(1,30)) 
    beta_true<-c(0.5,rep(-1,15),rep(0.5,15)) 
    
  } else if(i==16){
    b_true<-c(7,rep(1,30)) 
    beta_true<-c(0.5,rep(-1,15),rep(0.5,15)) 
    
  }
  
  mean.ppi<-rep(NA,nrep)
  mean.bag1<-rep(NA,nrep)
  
  for(j in 1:nrep){
    X_all<-BMIR2_genX(n=sim_param[i,1],m=sim_param[i,2],d=sim_param[i,3])
    X<-lapply(X_all, function(x){x[[1]]})
    delta_true<-lapply(X_all, function(x){x[[2]]})

    y = unlist(lapply(X_all, BMIR2_genY))
    
    prime_prop<-unlist(lapply(X_all, function(x){mean(x[[2]])}))
    mean.ppi[j]<-mean(prime_prop)
    mean.bag1[j]<-mean(y)
    
    save(X,y,beta_true,b_true, delta_true,file = paste0("simdat_",i,"_",j, ".Rdata"))
  }
  cat("Scenario:",i,"created. ",mean(mean.ppi),';',mean(mean.bag1),"\n")
}
