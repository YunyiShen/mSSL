// B_theta_update.h
// [[Rcpp::depends(RcppArmadillo)]]

// this is the helper function a naive implementation of cgLASSO with uniform L1 penalty, not particularly efficient
#ifndef CGLASSO_CGLASSO_B_UPDATE_H
#define CGLASSO_CGLASSO_B_UPDATE_H

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <RcppArmadillo.h>

#define EPS (double(2.22E-16))

using namespace Rcpp;
using namespace arma;
using namespace std;


namespace cgLASSO{



inline double cglassoobjective(int n, int p, int q,
                 const arma::mat &S,
                 const arma::mat &B,
                 const arma::mat &X,
                 const arma::mat &Omega,
                 const arma::mat &Sigma,
                 double lambda,
                 double xi,
                 int diag_penalty)
{

  double value;
  double sign;
  log_det(value, sign, Sigma);
  // this uses the data transforming trick. 
  arma::mat SOmega = S * Sigma;
  //arma::mat MSigma = B.t() * X.t() * X * B * Sigma / n;
  double log_like = value * sign - trace(SOmega);// - trace(MSigma);
  log_like *= n / 2;
  double log_pi_B = 0.0;
  double log_pi_Omega = 0.0;

  for (int j = 0; j < p; j++)
  {
    for (int k = 0; k < q; k++)
    {
      log_pi_B += (log(lambda) -1.0 * lambda * abs(B(j, k)));
    }
  }
  for (int k = 0; k < q; k++)
  {
    for (int kk = k + 1; kk < q; kk++)
    {
      log_pi_Omega += (log(xi) -1.0 * xi * abs(Omega(k, kk)));
    }
  }
  if (diag_penalty == 1)
  {
    for (int k = 0; k < q; k++)
    {
      log_pi_Omega += (log(xi) -1.0 * xi * abs(Omega(k, k)));
    }
  }


  double obj = log_like + log_pi_B + log_pi_Omega;

  return obj;
}


inline double lassoobjective(int n, int p, int q,arma::mat S, arma::mat B, arma::mat Omega, double lambda, double xi, int diag_penalty){
  
  double value;
  double sign;
  log_det(value, sign, Omega);
  
  arma::mat SOmega = S * Omega;
  double log_like = value*sign - trace(SOmega);
  log_like *= n/2;
  double log_pi_B = 0.0;
  double log_pi_Omega = 0.0;
  
  for(int j = 0; j < p; j++){
    for(int k = 0; k < q; k++){
      log_pi_B += (log(lambda) -1.0 * lambda * abs(B(j,k)));
    }
  }
  for(int k = 0; k < q; k++){
    for(int kk = k+1; kk < q; kk++){
      log_pi_Omega += (log(xi) - 1.0 * xi * abs(Omega(k,kk)));
    }
  }
  if(diag_penalty == 1){
    for(int k = 0; k < q; k++){
      log_pi_Omega += (log(xi) - 1.0 * xi * abs(Omega(k,k)));
    }
  }
  
  
  
  double obj = log_like + log_pi_B + log_pi_Omega ;
  
  
  return obj;
}



inline void B_coord_desc(const int n, const int p, const int q, arma::mat& B, arma::mat& R, arma::mat& tXR, arma::mat& S,
                  const arma::mat Omega,const arma::mat X, const arma::mat tXX,
                  const double lambda, const double xi,
                  const int max_iter, const double eps, const int verbose){
  // make a copy to be safe
  arma::mat B_init = B; 
  arma::mat R_init = R;
  arma::mat tXR_init = tXR;
  arma::mat S_init = S;
  double g = 0.0;
  double pstar0 = 0.0;
  double pstar0_inv = 0.0;
  double pstar = 0.0;
  double lambda_star = 0.0;
  double lambda_star0 = 0.0;
  double w_kk;
  double Delta = 0.0;
  double b_old = 0.0;
  double b_new = 0.0;
  double b_shift;
  double z = 0.0;
  double z_sgn = 1.0;
  int converged = 1;
  int violations = 0;
  int iter = 0;

  //double obj = 0.0;
  //double old_obj = 0.0;
  //int obj_counter = 0;
  
  arma::mat e1(p,q);
  e1.zeros(); // the active set
  for(int j = 0; j < p; j++){
    for(int k = 0; k < q; k++){
	  if(B(j,k) != 0) e1(j,k) = 1;
	  }
  }
  
    while(iter < max_iter){
      while(iter < max_iter){
        iter++;
        converged = 1;
        for(int j = 0; j < p; j++){
          for(int k = 0; k < q; k++){
            w_kk = Omega(k,k);
            if(e1(j,k) == 1){
              b_old = B(j,k);
              z = dot(tXR.row(j), Omega.row(k))/w_kk + n*b_old;
              z_sgn = (z > 0) - (z < 0);
              b_new = 1.0/n * z_sgn * max(abs(z) - lambda/w_kk, 0.0);
              B(j,k) = b_new;
              b_shift = b_old - b_new;
              R.col(k) += X.col(j) * b_shift;
              S.row(k) += b_shift * tXR.row(j)/n;
              S.col(k) += tXR.row(j).t() * b_shift/n;
              S(k,k) += tXX(j,j) * b_shift * b_shift/n;
              tXR.col(k) += tXX.col(j) * b_shift;
              if(abs(b_shift/b_old) > eps) converged = 0;
            }
          }
        }
        if(converged == 1) break;
      }
      violations = 0;
      for(int j = 0; j < p; j++){
        for(int k = 0; k < q; k++){
          w_kk = Omega(k,k);
          if(e1(j,k) == 0){
            b_old = B(j,k);
            z = dot(tXR.row(j), Omega.row(k))/w_kk + n*b_old;
            z_sgn = (z > 0) - (z < 0);
            b_new = 1.0/n * z_sgn * max(abs(z) - lambda/w_kk, 0.0);
            if(b_new != 0){ // violation found!
             violations++;
             B(j,k) = b_new;
             b_shift = b_old - b_new;
             R.col(k) += X.col(j) * b_shift;
             S.row(k) += b_shift * tXR.row(j)/n;
             S.col(k) += tXR.row(j).t() * b_shift/n;
             S(k,k) += tXX(j,j) * b_shift * b_shift/n;
             tXR.col(k) += tXX.col(j) * b_shift;
            }
          }
		      if(B(j,k) != 0) e1(j,k) = 1;
		      else e1(j,k) = 0;
        }
      }
      if(violations == 0) break;
    }
   
  if((iter == max_iter) & ((converged != 1 || violations != 0))){ // hit max iter and either: did not converge on active set or found violations in in-active set
    if(verbose == 1) Rcout << "        [B_coord_desc]: Did not converge. Violations = " << violations << endl;
  }
}



inline void update_B(const int n, const int p, const int q, arma::mat& B, arma::mat& R, arma::mat& tXR, arma::mat&S, const arma::mat Omega, 
                    const arma::mat X, const arma::mat tXX, const double lambda,  const double xi,  const int diag_penalty,
                    const int max_iter, const double eps, const int verbose){
  arma::mat B_old = B;
  arma::mat B_orig = B;
  arma::mat R_orig = R;
  arma::mat S_orig = S;
  arma::mat tXR_orig = tXR;
  //double theta_orig = theta;
  int b_max_iter = 5*max_iter;
  
  double obj = 0.0;
  double old_obj = 0.0;
  //int obj_counter = 0;
  
  int converged = 1;
  int iter = 0;
  while(iter < max_iter){
    iter++;
    B_old = B;
    old_obj = lassoobjective(n, p, q, S, B, Omega, lambda, xi, diag_penalty);
    B_coord_desc(n, p, q, B, R, tXR, S, Omega, X, tXX, lambda, xi,b_max_iter, eps,verbose);
    obj = lassoobjective(n,p,q,S, B, Omega, lambda, xi, diag_penalty);
    
    // check convergence
    converged = 1;
    
    for(int j = 0; j < p; j++){
      for(int k = 0; k < q; k++){
        if((B_old(j,k) == 0) & (B(j,k) != 0)) converged = 0;
        else if( (B_old(j,k) != 0) & (abs( (B_old(j,k) - B(j,k))/B_old(j,k)) > eps) ) converged = 0;
      }
    }
    
    if(verbose == 1) Rcout << "    [B_update]: iter = " << iter << " num B != 0: " << accu(B != 0) << endl;
    
    if(converged == 1) break;
  }
  if((iter == max_iter) & (converged == 0)){
    if(verbose == 1) Rcout << "Updating B failed to converge." << endl;
  }
  
}

}

#endif


