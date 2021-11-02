//--------------------------------------------------------------
// Header (header)
//--------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;


//--------------------------------------------------------------
// Functions (Functions_cpp)
//--------------------------------------------------------------
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}

double log_sum_exp(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  double maxv = x(maxi);
  if (!(maxv > -datum::inf)) {
    return -datum::inf;
  }
  double cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) & (x(i) > -datum::inf)) {
      cumsum += exp(x(i) - maxv);
    }
  }
  return maxv + log1p(cumsum);
}


double c_rtexp_onesided(double a)
{
  bool stop = false;
  double lambda = 0.5 * (a + sqrt(pow(a, 2) + 4));
  double z;
  
  while(!stop)
  {
    z = a - log(R::runif(0,1))/lambda;
    stop = (-2* log(R::runif(0,1)) > pow(z - lambda, 2));
  }
  return z;
}

// [[Rcpp::export]]
arma::vec ctruncnorm(arma::vec mean, 
                     arma::vec sd, 
                     bool upper){
  int m = mean.size();
  vec p = normcdf(- mean / sd);
  vec res = zeros(m);
  if(upper){
    // right truncated
    p += (1.0 - p) % randu(m);
    for(int j = 0; j < m; j++) {
      if(- mean(j) / sd(j) > 3.48672170399){
        res(j) = sd(j) * c_rtexp_onesided(- mean(j) / sd(j)) + mean(j);
      } else{
        res(j) = R::qnorm(p(j), mean(j), sd(j), true, false);
      }
    }
  } else{
    // left truncated
    p = p % randu(m);
    for(int j = 0; j < m; j++) {
      if(- mean(j) / sd(j) < -3.48672170399){
        res(j) = - sd(j) * c_rtexp_onesided(mean(j) / sd(j)) + mean(j);
      } else{
        res(j) = R::qnorm(p(j), mean(j), sd(j), true, false);
      }
    }
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::List MICProB_1Gibbs_cpp(
    arma::mat X,
    arma::vec y,
    arma::vec ninst,
    
    arma::vec hp_mu_beta,
    arma::vec hp_mu_b,
    arma::mat hp_Sig_beta,
    arma::mat hp_Sig_b,
    
    arma::vec beta,
    arma::vec b,
    arma::vec delta,
    arma::vec u,
    arma::vec z,
    
    arma::mat hp_Sig_beta_inv,
    arma::mat hp_Sig_b_inv,
    arma::mat V_b){
  
  int n = y.size(), d = X.n_cols;
  
  // set X_delta, mu_z
  mat X_delta = zeros(n, d);
  vec mu_z = zeros(n);
  int pos = 0;
  for(int i = 0; i < n; i++){
    if(accu(delta.subvec(pos, pos + ninst[i] - 1)) != 0){
      X_delta.row(i) = sum(X.rows(pos + find(delta.subvec(pos, pos + ninst[i] - 1) == 1)),
                  0); // column sum
    }
    pos += ninst[i];
  }
  mu_z = beta(0) + X_delta * beta.tail(d);
  
  // Rcout << "mu_z=\n" << mu_z << "\n";
  
  // update z
  colvec one_colvec = ones(n);
  //// y == 1
  uvec idx = find(y == 1);
  z(idx) = ctruncnorm(mu_z(idx), one_colvec(idx), true);
  //// y == 0
  idx = find(y == 0);
  z(idx) = ctruncnorm(mu_z(idx), one_colvec(idx), false);
  // z[y==1]<-rtruncnorm(sum(y==1),mean=mu_z[y==1],sd=1,a=0,b=Inf)
  //   z[y==0]<-rtruncnorm(sum(y==0),mean=mu_z[y==0],sd=1,a=-Inf,b=0)
  
  
  // update beta
  mat X1_delta = join_rows(one_colvec, X_delta);
  // X1_delta<-cbind(1,X_delta)
  mat V_beta = inv_sympd(hp_Sig_beta_inv + X1_delta.t() * X1_delta);
  vec m_beta = V_beta * (hp_Sig_beta_inv * hp_mu_beta + X1_delta.t() * z);
  beta = m_beta + (randn(1, d + 1) * chol(V_beta)).t();
  
  // update delta
  vec mu_u = b(0) + X * b.tail(d);
  vec probit_prob = normcdf(mu_u);
  
  // Rcout << "mu_u=\n" << mu_u << "\n";
  
  
  pos = 0;
  double A, B;
  for(int i = 0; i < n; i++){
    mat X_j = X.rows(pos, pos + ninst[i] - 1);
    vec delta_j = delta.subvec(pos, pos + ninst[i] - 1);
    for(int j = 0; j < ninst[i]; j++){
      double tmp = z(i) - beta(0) - 
        accu(X_j.rows(find(delta_j == 1)) * beta.tail(d)) +
        delta_j(j) * dot(X_j.row(j), beta.tail(d));
      
      A = exp(-0.5 * pow(tmp - dot(X_j.row(j), beta.tail(d)), 2.0));
      B = exp(-0.5 * pow(tmp, 2.0));
      
      double prim_prob = (A * probit_prob(pos + j)) / 
        (A * probit_prob(pos + j) + B * (1 - probit_prob(pos + j)));
      
      if(std::isnan(prim_prob)){
        prim_prob = 0;
      }
      // Rcout << "prim_prob=\n" << prim_prob << "\n";
      if(conv_to<double>::from(randu(1)) < prim_prob){
        delta(pos + j) = 1;
      } else{
        delta(pos + j) = 0;
      }  
    }
    pos += ninst[i];
  }
  
  // update u
  //// delta == 1
  vec one_colvec2 = ones(X.n_rows);
  uvec idx_delta1 = find(delta == 1);
  u(idx_delta1) = ctruncnorm(mu_u(idx_delta1), one_colvec2(idx_delta1), true);
  //// delta == 0
  uvec idx_delta2 = find(delta == 0);
  u(idx_delta2) = ctruncnorm(mu_u(idx_delta2), one_colvec2(idx_delta2), false);
  // u [delta==1]<-rtruncnorm(sum(delta==1),mean=mu_u[delta==1],sd=1,a=0,b=Inf);
  // u[delta==0]<-rtruncnorm(sum(delta==0),mean=mu_u[delta==0],sd=1,a=-Inf,b=0)

  // // update u
  // //// delta == 1
  // one_colvec = ones(X.n_rows);
  // uvec idx_delta1 = find(delta == 1);
  // u(idx) = ctruncnorm(mu_u(idx), one_colvec(idx), true);
  // //// delta == 0
  // idx = find(delta == 0);
  // u(idx) = ctruncnorm(mu_u(idx), one_colvec(idx), false);
  // // u [delta==1]<-rtruncnorm(sum(delta==1),mean=mu_u[delta==1],sd=1,a=0,b=Inf);
  // // u[delta==0]<-rtruncnorm(sum(delta==0),mean=mu_u[delta==0],sd=1,a=-Inf,b=0)
  
  // update b
  // one_colvec = ones(X.n_rows);
  vec m_b = V_b * (hp_Sig_b_inv * hp_mu_b + join_rows(one_colvec2, X).t() * u);
  b = m_b + (randn(1, d + 1) * chol(V_b)).t();
  // Rcout << "here3\n";
  
  return List::create(
    Named("beta") = beta,
    Named("b") = b,
    Named("delta") = delta,
    Named("u") = u,
    Named("z") = z
  );
}

// [[Rcpp::export]]
Rcpp::List Predict_cpp(
    arma::mat X_new,
    arma::vec ninst_new,
    arma::vec beta,
    arma::vec b){
  
  int n = ninst_new.size();
  int d = X_new.n_cols;
  // prediction for new bags
  vec prime_prob = normcdf(b(0) + X_new * b.tail(d));
  
  // generate 50 replicates of delta
  vec bag_prob = zeros(n);
  // mat delta_samps = zeros(ninst, n_rep_delta);
  
  vec u = randu(accu(ninst_new));
  vec delta = arma::conv_to<arma::vec>::from(prime_prob > u);
  mat X_new_delta = zeros(n, d);
  int pos = 0;
  for(int i = 0; i < n; i++){  
    int nprime = accu(delta.subvec(pos, pos + ninst_new[i] - 1));
    if(nprime != 0){
      X_new_delta.row(i) = sum(X_new.rows(pos + find(delta.subvec(pos, pos + ninst_new[i] - 1) == 1)),
                        0); // column sum
    }
    pos += ninst_new[i];
  }
  vec mu_z = beta(0) + X_new_delta * beta.tail(d);
  bag_prob += normcdf(mu_z);

  return List::create(
    Named("y_pred_prob") = bag_prob,
    Named("delta_pred_prob") = prime_prob
  );
  //   X_test_delta<-matrix(NA,n_test,d)
  //   mu_z<-rep(NA,n_test)
  //   bag_prob<-matrix(NA,n_test,dim(delta_samps)[2])
  //   
  //   for(k in 1:dim(delta_samps)[2]){
  //     for(i in 1:n_test){
  //       delta_test_i<-delta_samps[membership_test==i,k]
  //       X_test_i<-X_test[membership_test==i,,drop=F]
  //       
  //       if(sum(delta_test_i)==0){
  //         X_test_delta[i,]<-rep(0,d)
  //       } else{
  //         X_test_delta[i,]<-colSums(X_test_i[delta_test_i==1,,drop=F])
  //       }
  //       mu_z[i]<-beta[1] + X_test_delta[i,] %*% beta[-1]
  //       
  // #cat("k=",k," i=",i,"\n")
  //     }
  //     bag_prob[,k]<-pnorm(mu_z)
  //   }
  //   
  //   y_pred_prob[iter/nthin,]<-rowMeans(bag_prob)
  //     delta_pred_prob[iter/nthin,]<-prime_prob
}