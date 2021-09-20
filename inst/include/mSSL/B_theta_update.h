// B_theta_update.h
// [[Rcpp::depends(RcppArmadillo)]]


#ifndef MSSL_B_THETA_UPDATE_H
#define MSSL_B_THETA_UPDATE_H

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


namespace mSSL{

inline double gSSL_objective(int n, int q, arma::mat S, arma::mat Omega, double xi1, double xi0, double eta, int diag_penalty, arma::vec eta_hyper_params){
  
  double value;
  double sign;
  log_det(value, sign, Omega);
  
  arma::mat SOmega = S * Omega;
  double log_like = value*sign - trace(SOmega);
  log_like *= n/2;  
  double log_pi_Omega = 0.0;

  for(int k = 0; k < q; k++){
    for(int kk = k+1; kk < q; kk++){
      log_pi_Omega += log(eta * xi1 * exp(-1.0 * xi1 * abs(Omega(k,kk))) + (1.0 - eta) * xi0 * exp(-1.0 * xi0 * abs(Omega(k,kk))));
    }
  }
  if(diag_penalty == 1){
    for(int k = 0; k < q; k++){
      log_pi_Omega += log(xi1 * exp(-1.0 * xi1 * abs(Omega(k,k))));
    }
  }
  
  double a_eta = eta_hyper_params(0);
  double b_eta = eta_hyper_params(1);
  
  double log_pi_eta = (a_eta - 1.0) * log(eta) + (b_eta - 1.0) * log(1.0 - eta);
  
  double obj = log_like + log_pi_Omega + log_pi_eta;
  return obj;
}



inline double cgobjective(int n, int p, int q,
                 const arma::mat &S,
                 const arma::mat &B,
                 const arma::mat &X,
                 const arma::mat &Omega,
                 const arma::mat &Sigma,
                 double lambda1,
                 double lambda0,
                 double xi1,
                 double xi0, double theta, double eta,
                 int diag_penalty, arma::vec theta_hyper_params,
                 arma::vec eta_hyper_params)
{

  double value;
  double sign;
  log_det(value, sign, Omega);

  arma::mat SOmega = S * Omega;
  arma::mat MSigma = B.t() * X.t() * X * B * Sigma / n;
  double log_like = value * sign - trace(SOmega) - trace(MSigma);
  log_like *= n / 2;
  double log_pi_B = 0.0;
  double log_pi_Omega = 0.0;

  for (int j = 0; j < p; j++)
  {
    for (int k = 0; k < q; k++)
    {
      log_pi_B += log(theta * lambda1 * exp(-1.0 * lambda1 * abs(B(j, k))) + (1.0 - theta) * lambda0 * exp(-1.0 * lambda0 * abs(B(j, k))));
    }
  }
  for (int k = 0; k < q; k++)
  {
    for (int kk = k + 1; kk < q; kk++)
    {
      log_pi_Omega += log(eta * xi1 * exp(-1.0 * xi1 * abs(Omega(k, kk))) + (1.0 - eta) * xi0 * exp(-1.0 * xi0 * abs(Omega(k, kk))));
    }
  }
  if (diag_penalty == 1)
  {
    for (int k = 0; k < q; k++)
    {
      log_pi_Omega += log(xi1 * exp(-1.0 * xi1 * abs(Omega(k, k))));
    }
  }

  double a_theta = theta_hyper_params(0);
  double b_theta = theta_hyper_params(1);
  double a_eta = eta_hyper_params(0);
  double b_eta = eta_hyper_params(1);

  double log_pi_theta = (a_theta - 1.0) * log(theta) + (b_theta - 1.0) * log(1.0 - theta);
  double log_pi_eta = (a_eta - 1.0) * log(eta) + (b_eta - 1.0) * log(1.0 - eta);

  double obj = log_like + log_pi_B + log_pi_Omega + log_pi_theta + log_pi_eta;

  return obj;
}


// updated July 18: added indicator for diagonal penalty.
inline double objective(int n, int p, int q,arma::mat S, arma::mat B, arma::mat Omega, double lambda1, double lambda0, double xi1, double xi0, double theta, double eta, int diag_penalty, arma::vec theta_hyper_params, arma::vec eta_hyper_params){
  
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
      log_pi_B += log(theta * lambda1 * exp(-1.0 * lambda1 * abs(B(j,k))) + (1.0 - theta) * lambda0 * exp(-1.0 * lambda0 * abs(B(j,k))));
    }
  }
  for(int k = 0; k < q; k++){
    for(int kk = k+1; kk < q; kk++){
      log_pi_Omega += log(eta * xi1 * exp(-1.0 * xi1 * abs(Omega(k,kk))) + (1.0 - eta) * xi0 * exp(-1.0 * xi0 * abs(Omega(k,kk))));
    }
  }
  if(diag_penalty == 1){
    for(int k = 0; k < q; k++){
      log_pi_Omega += log(xi1 * exp(-1.0 * xi1 * abs(Omega(k,k))));
    }
  }
  
  double a_theta = theta_hyper_params(0);
  double b_theta = theta_hyper_params(1);
  double a_eta = eta_hyper_params(0);
  double b_eta = eta_hyper_params(1);
  
  double log_pi_theta = (a_theta - 1.0) * log(theta) + (b_theta - 1.0) * log(1.0 - theta);
  double log_pi_eta = (a_eta - 1.0) * log(eta) + (b_eta - 1.0) * log(1.0 - eta);
  
  double obj = log_like + log_pi_B + log_pi_Omega + log_pi_theta + log_pi_eta;
  
  
  return obj;
}



inline void B_coord_desc(const int n, const int p, const int q, arma::mat& B, arma::mat& R, arma::mat& tXR, arma::mat& S, const double theta,
                  const arma::mat Omega, const double eta,const arma::mat X, const arma::mat tXX,
                  const double lambda1, const double lambda0, const double xi1, const double xi0,
                  arma::vec theta_hyper_params, arma::vec eta_hyper_params, const int max_iter, const double eps, const int verbose){
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
  
  if(lambda1 == lambda0){
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
              b_new = 1.0/n * z_sgn * max(abs(z) - lambda1/w_kk, 0.0);
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
            b_new = 1.0/n * z_sgn * max(abs(z) - lambda1/w_kk, 0.0);
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
  } else{
    while(iter < max_iter){
      while(iter < max_iter){ // inner loop over the active set
        iter++;
        converged = 1;
        for(int j = 0; j < p; j++){
          for(int k = 0; k < q; k++){
            w_kk = Omega(k,k);
            if(e1(j,k) == 1){
              b_old = B(j,k);
              z = dot(tXR.row(j), Omega.row(k))/w_kk + n*b_old;
              z_sgn = (z > 0) - (z < 0);
              // get what we need for the Delta threshold

              pstar0_inv = 1.0 + (1.0 - theta)/theta * lambda0/lambda1;
              pstar0 = 1.0/pstar0_inv;
              
              lambda_star0 = lambda1 * pstar0 + lambda0 * (1.0 - pstar0);
              pstar = 1.0/(1.0 + (1.0 - theta)/theta * lambda0/lambda1 * exp(-1 * abs(b_old) * (lambda0 - lambda1)));
              lambda_star = lambda1 * pstar + lambda0 * (1.0 - pstar);
              // check whether g(0) > 0:
              g = (lambda_star0 - lambda1)*(lambda_star0 - lambda1) - 2.0 * n * w_kk * log(pstar0_inv);
              if((g > 0) && (lambda0 - lambda1 > sqrt(n)/(2.0 * sqrt(w_kk)))){ // now we can use the upper bound on Delta
                Delta = sqrt(2*n*log(pstar0_inv)/w_kk) + lambda1/w_kk;
              } else { // i.e. g(0) < 0, or g(0) > 0 but lambda1 - lambda0 not large enough. set Delta = lambda_star_0
                Delta = lambda_star0/w_kk;
              }          
              // now do refined thresholding
              if(abs(z) <= Delta){
                b_new = 0.0;
              } else{
                b_new = 1.0/n * z_sgn * max(abs(z) - lambda_star/w_kk, 0.0); 
              }

              b_shift = b_old - b_new;
              B(j,k) = b_new;
              R.col(k) += X.col(j) * b_shift;
              S.row(k) += b_shift * tXR.row(j)/n;
              S.col(k) += tXR.row(j).t() * b_shift/n;
              S(k,k) += tXX(j,j) * b_shift * b_shift/n;
              tXR.col(k) += tXX.col(j) * b_shift;
              if(abs(b_shift/b_old) > eps) converged = 0;
            }
          }
        } // finished looping over active set
        if(converged == 1) break;
      }// break out of inner loop after converging on active set
      violations = 0;
      for(int j = 0; j < p; j++){
        for(int k = 0; k < q; k++){
          w_kk = Omega(k,k);
          if(e1(j,k) == 0){ // check for violations on the in-active set
            b_old = B(j,k);
            z = dot(tXR.row(j), Omega.row(k))/w_kk + n*b_old;
            z_sgn = (z > 0) - (z < 0);
              // get what we need for the Delta threshold

            pstar0_inv = 1.0 + (1.0 - theta)/theta * lambda0/lambda1;
            pstar0 = 1.0/pstar0_inv;
              
            lambda_star0 = lambda1 * pstar0 + lambda0 * (1.0 - pstar0);
            pstar = 1.0/(1.0 + (1.0 - theta)/theta * lambda0/lambda1 * exp(-1 * abs(b_old) * (lambda0 - lambda1)));
            lambda_star = lambda1 * pstar + lambda0 * (1.0 - pstar);
              // check whether g(0) > 0:
            g = (lambda_star0 - lambda1)*(lambda_star0 - lambda1) - 2.0 * n * w_kk * log(pstar0_inv);
            if((g > 0) && (lambda0 - lambda1 > sqrt(n)/(2.0 * sqrt(w_kk)))){ // now we can use the upper bound on Delta
              Delta = sqrt(2*n*log(pstar0_inv)/w_kk) + lambda1/w_kk;
            } else { // i.e. g(0) < 0, or g(0) > 0 but lambda1 - lambda0 not large enough. set Delta = lambda_star_0
              Delta = lambda_star0/w_kk;
            }          
              // now do refined thresholding
            if(abs(z) <= Delta){
              b_new = 0.0;
            } else{
              b_new = 1.0/n * z_sgn * max(abs(z) - lambda_star/w_kk, 0.0); 
            }
            if(b_new != 0){ // found a violation!
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
          // update the active set at this time
          if(B(j,k) != 0) e1(j,k) = 1;
          else e1(j,k) = 0;
        }
      }
      if(violations == 0) break;
    }
  }
  if((iter == max_iter) & ((converged != 1 || violations != 0))){ // hit max iter and either: did not converge on active set or found violations in in-active set
    if(verbose == 1) Rcout << "        [B_coord_desc]: Did not converge. Violations = " << violations << endl;
  }
}


inline void update_theta(const int n, const int p, const int q, double& theta, const arma::mat B, const double lambda1, const double lambda0, arma::vec theta_hyper_params){
  double pstar = 0.0;
  double a_theta = theta_hyper_params(0);
  double b_theta = theta_hyper_params(1);
  
  double grad = 0.0;
  double hess = 0.0;
  double step = 1.0;
  int iter = 0;
  double theta_orig = theta;
  
  if(lambda0 == lambda1){
    theta = 1-1e-3;
  } else{
    while((iter < 1e6) & (abs(step) > 1e-3)){
      iter++;
      //Rcout << "iter = " << iter << " theta = " << theta << endl;
      grad = (a_theta - 1.0)/theta - (b_theta - 1.0)/(1.0 - theta);
      hess = -1.0 * (a_theta - 1.0)/(theta * theta) - (b_theta - 1.0)/((1.0 - theta)*(1.0 - theta));
      for(int j = 0; j < p; j++){
        for(int k = 0; k < q; k++){
          pstar = 1.0/(1.0 + (1.0 - theta)/theta * lambda0/lambda1 * exp(-1.0 * abs(B(j,k)) * (lambda0 - lambda1)));
          grad += pstar/theta - (1.0 - pstar)/(1.0 - theta);
          hess -= (pstar/theta - (1.0 - pstar)/(1.0 - theta))*(pstar/theta - (1.0 - pstar)/(1.0 - theta));
        }
      }
      step = grad/hess;
    //Rcout << "step = " << step << "theta = " << theta << endl;
      if(step > theta){ // full Newton step would take us outside the interval [0,1]
        step = theta - 1e-3; // this sets next value of theta to be 1e-3
      //Rcout << "new step = " << step << endl;
      } else if(step < theta - 1){ // full Newton step takes us outside the interval [0,1]
        step = theta + 1e-3 - 1; // this sets next value of theta to be 1 - 1e-6
      //Rcout << "new step =" << step << endl;
      }
      theta -= step;
    }
    if( (iter == 1e6) & (abs(step) > 1e-3)){
      // have not yet converged
      // eventually we should throw an error
    //  Rcout << "theta did not converge!" << endl;
      theta = theta_orig; // just set it back to the original value.
    }
  }
}

inline void update_B_theta(const int n, const int p, const int q, arma::mat& B, arma::mat& R, arma::mat& tXR, arma::mat&S, double& theta, const arma::mat Omega, const double eta, 
                    const arma::mat X, const arma::mat tXX, const double lambda1, const double lambda0, const double xi1, const double xi0, const int diag_penalty,
                    arma::vec theta_hyper_params, arma::vec eta_hyper_params, const int max_iter, const double eps, const int verbose){
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
    old_obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty,theta_hyper_params, eta_hyper_params);
    B_coord_desc(n, p, q, B, R, tXR, S, theta, Omega, eta, X, tXX, lambda1, lambda0, xi1, xi0,theta_hyper_params, eta_hyper_params,b_max_iter, eps,verbose);
    update_theta(n, p, q, theta, B, lambda1, lambda0, theta_hyper_params);
    obj = objective(n,p,q,S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty,theta_hyper_params, eta_hyper_params);
    
    // check convergence
    converged = 1;
    
    for(int j = 0; j < p; j++){
      for(int k = 0; k < q; k++){
        if((B_old(j,k) == 0) & (B(j,k) != 0)) converged = 0;
        else if( (B_old(j,k) != 0) & (abs( (B_old(j,k) - B(j,k))/B_old(j,k)) > eps) ) converged = 0;
      }
    }
    
    if(verbose == 1) Rcout << "    [B_theta_update]: iter = " << iter << " num B != 0: " << accu(B != 0) << endl;
    
//    if( (obj - old_obj)/old_obj < 1e-3) obj_counter++;
//    else obj_counter = 0;
//    if((obj_counter == 5) & (converged == 0)){
//      Rcout << "[update_B_theta]: obj_counter == 5 & converged = 0" << endl;
//      converged = 1;
//    }
    if(converged == 1) break;
  }
  if((iter == max_iter) & (converged == 0)){
    if(verbose == 1) Rcout << "Updating B and theta failed to converge." << endl;
    //B = B_orig;
    //R = R_orig;
    //S = S_orig;
    //tXR = tXR_orig;
    //theta = theta_orig;
  }
  
}

}

#endif


