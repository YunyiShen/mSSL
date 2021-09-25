#ifndef MSSL_MDCPE_H
#define MSSL_MDCPE_H

#include <RcppArmadillo.h>
#include <linconGaussR.h>
#include <QUIC.h>
#include <workingparam.h>
#include "B_theta_update.h"

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace workingparam;
using namespace linconGaussR;
using namespace workingparam;


namespace mSSL{

template<class workpara>
List mSSL_dcpe(workpara instance,
               List lambdas,
               List xis,
               arma::vec theta_hyper_params,
               arma::vec eta_hyper_params,
               int diag_penalty,
               int max_iter,
               double eps,
               int verbose,
               int n_rep = 1000,
              int nskp = 5)
{

  int n = instance.Y.n_rows;
  int q = instance.Y.n_cols;
  int p = instance.instance.X.n_cols;
  
  // Center columns of instance.X and Y
  // re-scale the re-centered columns of instance.X to have norm sqrt(n)
  
  double tmp_mu_x = 0.0;
  double tmp_weight_x = 0.0;
  arma::vec x_col_weights(p);
  arma::vec mu_x(p);
  for(int j = 0; j < p ; j++){
    tmp_mu_x = mean(instance.X.col(j));
    instance.X.col(j) -= tmp_mu_x;
    tmp_weight_x = norm(instance.X.col(j))/sqrt(n);
    instance.X.col(j) /= tmp_weight_x;
    mu_x(j) = tmp_mu_x;
    x_col_weights(j) = tmp_weight_x;
  }
  arma::vec mu_old(q,fill::zeros);
  int quic_max_iter = 5*max_iter;

  /*
  int max_iter = control_params["max_iter"];
  double eps = control_params["eps"];
  */
  double lambda1 = lambdas["lambda1"];
  arma::vec lambda_spike = lambdas["lambda0"];
  double lambda0 = lambda_spike(0);
  int L = lambda_spike.n_elem;
  
  double xi1 = xis["xi1"];
  arma::vec xi_spike = xis["xi0"];
  double xi0 = xi_spike(0);
  
  double eta = 0.5;
  double theta = 0.5;
  
  // initialize our parameters
  arma::vec alpha(q); // will hold the intercept. only computed at the end of all of the loops
  
  arma::mat B(p,q);
  B.zeros();
  arma::mat B_old = B; // used internally
  arma::mat tmp_B = B; // we will re-scale the working B matrix when we save the path
  arma::mat Omega(q,q);
  Omega.eye();
  arma::mat Omega_old = Omega; // used internally to check convergence in the EM algorithm
  arma::mat Sigma(q,q);
  Sigma.eye();
  arma::mat q_star(q,q);
  q_star.fill(0.5);
  arma::mat xi_star(q,q);
  xi_star.fill(xi1);
  // June 27: penalize diagonals
  xi_star.diag().fill(xi1);
  xi_star /= n; // the re-scaling is for compatability with QUIC.
  //xi_star.diag().zeros(); // for now, don't penalize the diagonals
  double tmp; //holds value of q_star
  int omega_non_zero = 0; // counter for the number of non-zero off-diagonal elements of Omega's
  
  if(verbose == 1) Rcout << "Initialized R, tXR, S" << endl;
  
  
  // if we ever fail to converge, we should re-set the values of B and Omega
  arma::mat B_reset = B;
  arma::mat Omega_reset = Omega;
  arma::mat Sigma_reset = Sigma;
  workpara instance_reset = instance;
  double theta_reset = theta;
  double eta_reset = eta;
  
  //Rcout << "preparing QUIC working parameters" << endl;
  // initialize stuff for QUIC
  cube res_quic(q,q,2); // cube to hold values from QUIC. 1st one is for Omega, 2nd is for Sigma
  
  double a_eta = eta_hyper_params(0);
  double b_eta = eta_hyper_params(1);
  
  
  int converged = 1;
  int iter = 0;
  
  // Store the conditional modes
  cube B0_path(p,q,L);
  cube Omega0_path(q,q,L);
  arma::vec theta0_path(L);
  arma::vec eta0_path(L);
  if(verbose == 1) Rcout << "Starting Step 1: Estimating conditional modes of B" << endl;
  time_t tp;
  int time_start = time(&tp);
  for(int l = 0; l < L; l++){
    lambda0 = lambda_spike(l);
    if(verbose == 1) Rcout << setprecision(6) << "Starting lambda0 = " << lambda0 << " num B non-zero = " << arma::accu(B != 0) << " theta = " << theta << endl;
    instance.update(mu_old,B,Sigma,n_rep,nskp);
    mu_old = instance.mu;
    update_B_theta(instance.n_B,p,q,B,instance.R,instance.tXR,instance.S,theta,Omega,eta,instance.X,instance.tXX,lambda1,lambda0,xi1,xi0,diag_penalty,theta_hyper_params,eta_hyper_params,max_iter,eps, verbose);
    tmp_B = B;
    tmp_B.each_col()/x_col_weights;
    B0_path.slice(l) = tmp_B;
    theta0_path(l) = theta;
  }
  if(verbose == 1) Rcout << "starting step 2 now" << endl;
  for(int l = 0; l < L; l++){
    xi0 = xi_spike(l);
    omega_non_zero = 0;
    for(int k = 0; k < q; k++){
      for(int kk = k+1; kk < q; kk++){
        if(Omega(k,kk) != 0) omega_non_zero++;
      }
    }
    if(verbose == 1) Rcout << "Starting xi0 = " << xi0 << " num Omega non-zero = " << omega_non_zero << " eta = " << eta << endl;
    iter = 0;
    while(iter < max_iter){
      iter++;
      Omega_old = Omega;
      // E-step
      if(xi0 == xi1){
        q_star.fill(1.0);
      } else{
        for(int k = 0; k < q; k++){
          for(int kk = k+1; kk < q; kk++){
            tmp = 1.0/(1.0 + (1.0-eta)/eta * xi0/xi1 * exp(-1.0 * abs(Omega(k,kk)) * (xi0 - xi1)));
            q_star(k,kk) = tmp;
            q_star(kk,k) = tmp;
          }
        }
      }
      q_star.diag().zeros();
      xi_star = xi1 * q_star + xi0 * (1 - q_star);
      if(diag_penalty == 1) xi_star.diag().fill(xi1);
      else if(diag_penalty == 0) xi_star.diag().fill(0);
      // M-step update of eta
      eta = (a_eta - 1 + arma::accu(q_star)/2)/(a_eta + b_eta -2 + q*(q-1)/2);
      // M-step update of Omega
      xi_star /= n; // QUIC needs everything to be scaled
      res_quic = my_quic(q, instance.S, xi_star, eps, quic_max_iter);
      Omega = res_quic.slice(0);
      Sigma = res_quic.slice(1);
      instance.postprocessing(B,Sigma,Omega);

      // check convergence
      converged = 1;
      omega_non_zero = 0;
      for(int k = 0; k < q; k++){
        for(int kk = k+1; kk<q; kk++){
          if( (Omega_old(k,kk) == 0) & (Omega(k,kk) != 0)) converged = 0;
          else if( (Omega_old(k,kk) != 0) & (abs((Omega_old(k,kk) - Omega(k,kk))/Omega_old(k,kk)) > eps))converged = 0;
          if(Omega(k,kk) != 0) omega_non_zero ++;
        }
      }
      if(converged == 1) break;
    } // closes while loop (i.e. the main EM)
    if( (iter == max_iter) & (converged == 0) ){
      Rcout << "Omega did not converge. Re-setting the values" << endl;
      Omega = Omega_reset;
      Sigma = Sigma_reset;
      eta = eta_reset;
    } else if(converged == 1){
      Omega_reset = Omega;
      Sigma_reset = Sigma;
      eta_reset = eta;
    }
    Omega0_path.slice(l) = Omega;
    
  } // closes for loop over l
  
  if(verbose == 1) Rcout << "starting Step 3" << endl;
  B_reset = B;
  Omega_reset = Omega;
  Sigma_reset = Sigma;
  instance_reset = instance;
  theta_reset = theta;
  eta_reset = eta;
  
  lambda0 = lambda_spike(L-1);
  xi0 = xi_spike(L-1);
  iter = 0;
  converged = 1;
  while(iter < max_iter){
    iter++;
    B_old = B;
    Omega_old = Omega;
    
    instance.update(mu_old, B,Sigma,n_rep,nskp);
    mu_old = instance.mu;
    update_B_theta(instance.n_B,p,q,B,instance.R,instance.tXR,instance.S,theta,Omega,eta,instance.X,instance.tXX,lambda1,lambda0,xi1,xi0,diag_penalty,theta_hyper_params,eta_hyper_params,max_iter,eps, verbose);
    
    // now update Omega and eta
    for(int k = 0; k < q; k++){
      for(int kk = k+1; kk < q; kk++){
        tmp = 1.0/(1.0 + (1.0-eta)/eta * xi0/xi1 * exp(-1.0 * abs(Omega(k,kk)) * (xi0 - xi1)));
        q_star(k,kk) = tmp;
        q_star(kk,k) = tmp;
      }
    }
    q_star.diag().zeros();
    xi_star = xi1 * q_star + xi0 * (1 - q_star);
    if(diag_penalty == 1) xi_star.diag().fill(xi1);
    else if(diag_penalty == 0) xi_star.diag().fill(0);
    // M-step update of eta
    eta = (a_eta - 1 + arma::accu(q_star)/2)/(a_eta + b_eta -2 + q*(q-1)/2);
    // M-step update of Omega
    xi_star /= n;
    res_quic = my_quic(q, instance.S, xi_star, eps, quic_max_iter);
    Omega = res_quic.slice(0);
    Sigma = res_quic.slice(1);
    instance.postprocessing(B,Sigma,Omega);
    
    converged = 1;
    for(int j = 0; j < p; j++){
      for(int k = 0; k < q; k++){
        if( (B_old(j,k) == 0) & (B(j,k) != 0)) converged = 0;
        else if( (B_old(j,k) != 0) & (abs((B_old(j,k) - B(j,k))/B_old(j,k)) > eps)) converged = 0;
      }
    }
    omega_non_zero = 0;
    for(int k = 0; k < q; k++){
      for(int kk = k+1; kk < q; kk++){
        if( (Omega_old(k,kk) == 0) & (Omega(k,kk) != 0)) converged = 0;
        else if( (Omega_old(k,kk) != 0) & (abs( (Omega_old(k,kk) - Omega(k,kk))/Omega_old(k,kk)) > eps)) converged = 0;
        if(Omega(k,kk) != 0) omega_non_zero++;
      }
    }
    if(verbose == 1) Rcout << "iter = " << iter << " num nonzero B: " << arma::accu(B!=0) << " num nonzero Omega: " << omega_non_zero << " theta = " << theta << " eta = " << eta << endl;
    if(converged == 1) break;
  }
  if( (iter == max_iter) & (converged == 0) ){
    // did not converge!
    Rcout << "Omega and B did not converge" << endl;
  }
  
  int time_end = time(&tp);
  
  List results;
  tmp_B = B;
  tmp_B.each_col() /= x_col_weights;
  alpha = instance.mu - tmp_B.t() * mu_x;
  results["alpha"] = alpha;
  results["B"] = tmp_B;
  results["Omega"] = Omega;
  //results["obj"] = objective(n, p, q, S, B, Omega, lambda1, lambda_spike(L-1), xi1, xi_spike(L-1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
  results["B0_path"] = B0_path;
  //results["theta0_path"] = theta0_path;
  results["Omega0_path"] = Omega0_path;
  results["theta"] = theta;
  results["eta"] = eta;
  results["time"] = time_end - time_start;
  return results;
}

}


#endif