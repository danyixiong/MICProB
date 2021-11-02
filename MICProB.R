# library(coda)

#### Get hyperparameters ####
getHyperPars <- function(tidydata){
  D = tidydata$nfeature_inst
  res = list(
    #hp_mu_beta=c(mean(tidydata$label),rep(0,D)),
    #hp_mu_b=c(mean(tidydata$label),rep(0,D)),

    hp_mu_beta=rep(0, D + 1),
    hp_mu_b=rep(0,D + 1),
    
    #hp_Sig_beta=diag(D+1)/(D+1),
    #hp_Sig_b=diag(D+1)/(D+1)
    
    hp_Sig_beta = diag(10,D+1), 
    hp_Sig_b = diag(10,D+1)
    
    # for real data
    #hp_Sig_beta = diag(c(16,rep(4,D)),D+1), 
    #hp_Sig_b = diag(c(16,rep(4,D)),D+1)
  )
  return(res)
}

#### Get initial values ####
getInits <- function(tidydata,hyperpars){
  #beta = unlist(lapply(hyperpars$hp_mu_beta, function(mu_beta) rnorm(1, mu_beta, 10)))
  #b = unlist(lapply(hyperpars$hp_mu_b, function(mu_b) rnorm(1, mu_b, 10)))
  #delta = unlist(lapply(1:tidydata$nsample,function(i){rbinom(tidydata$ninst[i],1,mean(tidydata$label))}))
  
  ## use true parameter values
  beta = beta_true
  b = b_true
  delta = unlist(delta_true)
  #cat("Using true initials!\n")
  
  res = list(
    beta = beta,
    b = b,
    delta = delta
  )
  return(res)
}


#### Summarize all input data and parameters for mcmc chain #### 
getInputPars <- function(tidydata){
  
  list_hyperpars <- getHyperPars(tidydata)
  list_inits <- getInits(tidydata,list_hyperpars)
  
  res = list(
    ## data
    n = tidydata$nsample, # number of bags
    d = tidydata$nfeature_inst, # number of features
    m = tidydata$ninst, # number of instances per bag
    membership = tidydata$membership, # membership for instances
    y = tidydata$label, # bag labels
    X1 = model.matrix(~1 + ., data = Reduce(rbind, tidydata$feature_inst)), # design matrix
    ## hyperparameters 
    hp_mu_beta = list_hyperpars$hp_mu_beta,
    hp_mu_b = list_hyperpars$hp_mu_b,
    hp_Sig_beta = list_hyperpars$hp_Sig_beta,
    hp_Sig_b = list_hyperpars$hp_Sig_b,
    ## model parameters
    beta = list_inits$beta,
    b = list_inits$b,
    delta = list_inits$delta
  )
  return(res)
}


#### Fitting BMIR2 model ####

#### 1 Gibbs iteration in Rcpp ####
MICProB_sampler<-function(tidytrain,
                        tidytest,
                        ntotal,
                        nwarm,
                        nthin,
                        nchain,
                        #scale,
                        return_delta){
  # tidytrain=tidy_train
  # tidytest=tidy_train
  # ntotal=1000
  # nwarm=500
  # nthin=10
  # scale=F
  # return_delta=T
  
  cat("=============================================================\n")
  cat(sprintf("Probit Bayesian Multiple Instance Classification\n"))
  
  res_mcmc <- vector("list", nchain)
  
  for(nc in 1:nchain){
    
    # begin time
    start_time <- Sys.time()
    
    parlist <- getInputPars(tidytrain)
    
    # whether to scale the covariate matrix
    # if(scale){
    #   X1<-cbind(1,scale(parlist$X1[,-1]))
    #   #X1<-cbind(1,scale(parlist$X1[,-c(1,7)]), parlist$X1[,7]) # mute_type is categorical
    # } else{
    #   X1<-parlist$X1
    # }
    
    y<-parlist$y
    n<-parlist$n
    d<-parlist$d
    m<-parlist$m
    N<-sum(m)
    membership<-parlist$membership
    
    hp_mu_beta<-parlist$hp_mu_beta
    hp_mu_b<-parlist$hp_mu_b
    hp_Sig_beta<-parlist$hp_Sig_beta
    hp_Sig_b<-parlist$hp_Sig_b
    
    beta<-parlist$beta
    b<-parlist$b
    delta<-parlist$delta
    u = rep(0,N)
    z = rep(0,n)
    
    hp_Sig_beta_inv<-solve(hp_Sig_beta)
    hp_Sig_b_inv<-solve(hp_Sig_b)
    
    # posterior variance of b
    X1 <- parlist$X1
    V_b <- solve(hp_Sig_b_inv + crossprod(X1, X1))
    
    #cat("=============================================================\n")
    #cat(sprintf("Bayesian Multiple Instance Regression: chain" ,nc, " \n"))
    
    # begin time
    # start_time <- Sys.time()
    
    tick = 0.2
    
    # Gibbs sampling (warming up)
    
    cat("=============================================================\n")
    cat("Start warming up",nwarm,"MCMC samples!\n")
    cat("Progress: ")
    
    for(iter in 1:nwarm){
      if(iter %in% seq(round(tick*nwarm),nwarm,by=round(tick*nwarm))){
        cat(100*iter/nwarm,"% ...")
      }
      mcmc_res <- MICProB_1Gibbs_cpp(X = X1[,-1],y = y,
                                      ninst = m,
                                      hp_mu_beta = hp_mu_beta,
                                      hp_mu_b,
                                      hp_Sig_beta,
                                      hp_Sig_b,
                                      beta,
                                      b,
                                      delta,
                                      u,
                                      z,
                                      hp_Sig_beta_inv,
                                      hp_Sig_b_inv,
                                      V_b)
      # mcmc_res<-BMIR2_1Gibbs(X1,y,n,d,m,N,membership
      #                        ,hp_mu_beta,hp_mu_b,hp_Sig_beta,hp_Sig_b
      #                        ,beta,b,delta,u,z
      #                        ,hp_Sig_beta_inv,hp_Sig_b_inv,V_b)
      # update parameters
      beta = mcmc_res$beta
      b = mcmc_res$b
      delta = mcmc_res$delta
      u = mcmc_res$u
      z = mcmc_res$z
      
    } # end warm-up
    cat("\n")
    cat("Finish warming up!\n")
    cat("-------------------------------------------------------------\n")
    
    niter = ntotal - nwarm
    nsave = 1 + floor((niter - 1) /nthin)
    
    # posterior quantities to be saved
    beta_post<-matrix(NA,nrow=nsave,ncol=length(beta))
    b_post<-matrix(NA,nrow=nsave,ncol=length(b))
    delta_post<-matrix(NA,nrow=nsave,ncol=length(delta))
    
    pip_1chain<-rep(0,length(delta))
    mcmc_1chain <- list()
    
    # information from new bags
    X_test<-as.matrix(do.call(rbind,tidytest$feature_inst))
    ninst_test = tidytest$ninst
    #n_test<-tidytest$nsample
    #membership_test<-tidytest$membership
    
    # y_pred_prob and delta_pred_prob
    y_pred_prob<-matrix(NA,nrow=nsave,ncol=length(ninst_test))
    delta_pred_prob<-matrix(NA,nrow=nsave,ncol=nrow(X_test))
    
    cat("Start extracting",niter,"MCMC samples!\n")
    cat("Progress :")
    for(iter in 1:niter){
      if(iter %in% seq(round(tick*niter),niter,by=round(tick*niter))){
        cat(100*iter/niter,"% ...")
      }
      mcmc_res <- MICProB_1Gibbs_cpp(X = X1[,-1],y = y,
                                      ninst = m,
                                      hp_mu_beta = hp_mu_beta,
                                      hp_mu_b,
                                      hp_Sig_beta,
                                      hp_Sig_b,
                                      beta,
                                      b,
                                      delta,
                                      u,
                                      z,
                                      hp_Sig_beta_inv,
                                      hp_Sig_b_inv,
                                      V_b)
      # mcmc_res<-BMIR2_1Gibbs(X1,y,n,d,m,N,membership
      #                        ,hp_mu_beta,hp_mu_b,hp_Sig_beta,hp_Sig_b
      #                        ,beta,b,delta,u,z
      #                        ,hp_Sig_beta_inv,hp_Sig_b_inv,V_b)
      # update parmaeters
      beta = mcmc_res$beta
      b = mcmc_res$b
      delta = mcmc_res$delta
      u = mcmc_res$u
      z = mcmc_res$z
    
      # prediction for new bags
      #pred_res = Predict_cpp(X_test,ninst_test,beta,b)
      
      # prime_prob<-pnorm(b[1] + X_test%*%b[-1])
      # delta_samps<-replicate(50,rbinom(nrow(X_test),1,prime_prob))
      # 
      # X_test_delta<-matrix(NA,n_test,d)
      # mu_z<-rep(NA,n_test)
      # bag_prob<-matrix(NA,n_test,dim(delta_samps)[2])
      # 
      # for(k in 1:dim(delta_samps)[2]){
      #   for(i in 1:n_test){
      #     delta_test_i<-delta_samps[membership_test==i,k]
      #     X_test_i<-X_test[membership_test==i,,drop=F]
      #     
      #     if(sum(delta_test_i)==0){
      #       X_test_delta[i,]<-rep(0,d)
      #     } else{
      #       # sum model
      #       X_test_delta[i,]<-colSums(X_test_i[delta_test_i==1,,drop=F])
      #       # average model
      #       #X_test_delta[i,]<-colSums(X_test_i[delta_test_i==1,,drop=F])/sum(delta_test_i)
      #     }
      #     mu_z[i]<-beta[1] + X_test_delta[i,] %*% beta[-1]
      #     
      #     #cat("k=",k," i=",i,"\n")
      #   }
      #   bag_prob[,k]<-pnorm(mu_z)
      # }
      
      #y_pred_prob[iter,]<-pred_res$y_pred_prob
      #delta_pred_prob[iter,]<-pred_res$delta_pred_prob
      
      # save posterior samples and predicted bag probs
      if(iter %in% seq(nthin,niter,by=nthin)){ # thinning delta
        if(return_delta){
          delta_post[iter/nthin,]<-delta
        }
        pip_1chain = pip_1chain + delta
        beta_post[iter/nthin,]<-beta
        b_post[iter/nthin,]<-b
        
        pred_res = Predict_cpp(X_test,ninst_test,beta,b)
        y_pred_prob[iter/nthin,]<-pred_res$y_pred_prob
        delta_pred_prob[iter/nthin,]<-pred_res$delta_pred_prob
      
      }
      
    } # end extracting posterior samples
    
    pip_1chain = pip_1chain / nsave
    cat("\n")
    cat("Finish MCMC sampling!\n")
    cat("=============================================================\n")
    
    # elapsed time
    cat(sprintf("Elapsed time for chain%d=%.3f mins: MCMC sampling is done!\n", nc, difftime(Sys.time(), start_time, units = "mins")))
    
    # output
    mcmc_1chain[["beta"]]<-rbind(parlist$beta,beta_post)
    mcmc_1chain[["b"]]<-rbind(parlist$b,b_post)
    mcmc_1chain[["pip"]]<-pip_1chain
    
    if(return_delta){
      mcmc_1chain[["delta"]]<-rbind(parlist$delta, delta_post)
    } else{
      mcmc_1chain[["delta"]]<-NULL
    }
  
    mcmc_1chain[["y_pred"]]<-colMeans(y_pred_prob)
    mcmc_1chain[["delta_pred"]]<-colMeans(delta_pred_prob)
    
    res_mcmc[[nc]]<-mcmc_1chain
    
  }
  
  return(res_mcmc)
}
