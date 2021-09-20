//
//  mSSL_new.cpp
//  Source code for multivariate spike-and-slab (mSSL)
//  Includes function for both dynamic posterior exploration (DPE) and
// Dynamic conditonal posterior exploration (DCPE), as well as an
// EM algorithm for fitting Gaussian graphical models with SSLasso
//  Created by Sameer Deshpande on 1/25/19.
//


#include <RcppArmadillo.h>
#include <math.h>
#include <stdio.h>
#include <ctime>
#include <mSSL.h>
#include <QUIC.h>

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace mSSL;
using namespace quic;

// [[Rcpp::export]]
List mSSL_dpe(arma::mat X,
              arma::mat Y,
              List lambdas,
              List xis,
              arma::vec theta_hyper_params,
              arma::vec eta_hyper_params,
              int diag_penalty,
              int max_iter,
              double eps,
              int s_max_condition,
              int obj_counter_max,
              int verbose)
{
  int n = Y.n_rows;
  int q = Y.n_cols;
  int p = X.n_cols;
  
  // Center columns of X and Y
  // re-scale the re-centered columns of X to have norm sqrt(n)
  
  double tmp_mu_x = 0.0;
  double tmp_weight_x = 0.0;
  double tmp_mu_y = 0.0;
  arma::vec x_col_weights(p);
  arma::vec mu_x(p);
  arma::vec mu_y(q);
  for(int j = 0; j < p ; j++){
    tmp_mu_x = mean(X.col(j));
    X.col(j) -= tmp_mu_x;
    tmp_weight_x = norm(X.col(j))/sqrt(n);
    X.col(j) /= tmp_weight_x;
    mu_x(j) = tmp_mu_x;
    x_col_weights(j) = tmp_weight_x;
  }
  for(int k = 0; k < q; k++){
    tmp_mu_y = mean(Y.col(k));
    Y.col(k) -= tmp_mu_y;
    mu_y(k) = tmp_mu_y;
  }
  int quic_max_iter = 5*max_iter; // allow QUIC to run for more iterations

  /*
  // Read the control parameters
  int max_iter = control_params["max_iter"];
  double eps = control_params["eps"]; // typically 1e-3 but can be changed
  double s_max_condition = control_params["s_max_condition"]; // by default it should be 10*n
  int obj_counter_max = control_params["obj_counter_max"]; // maximum number of times that objective increases by < eps before terminating early
  */
  // Initialize grid of lambda and xi values
  double lambda1 = lambdas["lambda1"];
  arma::vec lambda_spike = lambdas["lambda0"];
  double lambda0 = lambda_spike(0);
  int L = lambda_spike.n_elem;
  double xi1 = xis["xi1"];
  arma::vec xi_spike = xis["xi0"];
  double xi0 = xi_spike(0);
  
  // Read hyper-parameters
  double a_eta = eta_hyper_params(0);
  double b_eta = eta_hyper_params(1);
  
  // Initialize the working parameters
  arma::vec alpha(q); // will hold the intercept. only computed at the end of all of the loops
  arma::mat B(p,q);
  B.zeros();
  arma::mat B_old = B; // used internally to check convergence
  arma::mat tmp_B = B; // used at the end to re-scale columns
  arma::mat Omega(q,q);
  Omega.eye();
  arma::mat Omega_old = Omega;
  arma::mat Sigma(q,q);
  Sigma.eye();
  double theta = 0.5;
  double eta = 0.5;
  double obj = 0.0; // value of objective function
  double old_obj = 0.0; // used to check convergence
  
  // Initialize the residual matrices and S
  arma::mat tYY = Y.t() * Y;
  arma::mat tXX = X.t() * X;
  arma::mat tXY = X.t() * Y;
  arma::mat R = Y - X * B;
  arma::mat tXR = X.t() * R;
  arma::mat tRR = R.t() * R;
  arma::mat S = tRR/n;
  arma::vec s_eval = eig_sym(S); // vector of eigenvalues of S
  
  // Initialize a copy for the re-starts
  arma::mat B_restart = B;
  arma::mat Omega_restart = Omega;
  arma::mat R_restart = R;
  arma::mat tXR_restart = tXR;
  arma::mat tRR_restart = tRR;
  arma::mat S_restart = S;
  double theta_restart = theta;
  double eta_restart = eta;
  
  // Initialize p_star and q_star and penalty matrices
  arma::mat q_star(q,q);
  q_star.fill(0.5);
  arma::mat xi_star(q,q);
  xi_star.fill(xi1/n);
  xi_star.diag().zeros();
  double tmp = 0.0; // holds the value of q_star(k,k')
  int omega_non_zero = 0; // counter for number of non-zero off-diagonal entries in Omega
  
  // Initialize stuff for QUIC
  cube res_quic(q,q,2); // cube to hold value of Omega and Sigma computed in QUIC
  
  
  
  // We need copies of B, Omega, R, tXR, S, theta, and eta that are propogated along the grid
  arma::mat B_left = B;
  arma::mat Omega_left = Omega;
  arma::mat Sigma_left = Sigma;
  arma::mat R_left = R;
  arma::mat tXR_left = tXR;
  arma::mat S_left = S;
  double theta_left = theta;
  double eta_left = eta;
  int left_index = 0;
  double obj_left = 0.0;
  
  arma::mat B_down = B;
  arma::mat Omega_down = Omega;
  arma::mat Sigma_down = Sigma;
  arma::mat R_down = R;
  arma::mat tXR_down = tXR;
  arma::mat S_down = S;
  double theta_down = theta;
  double eta_down = eta;
  int down_index = 0;
  double obj_down = 0.0;
  
  arma::mat B_diag = B;
  arma::mat Omega_diag = Omega;
  arma::mat Sigma_diag = Sigma;
  arma::mat R_diag = R;
  arma::mat tXR_diag = tXR;
  arma::mat S_diag = S;
  double theta_diag = theta;
  double eta_diag = eta;
  int diag_index = 0;
  double obj_diag = 0.0;
  
  // Initialize iterators, counters, etc.
  int s = 0;
  int t = 0;
  int current_index = 0;
  int obj_counter = 0; // counts number of times objective changes by less than 100*eps%
  int iter = 0;
  int converged = 0;
  int early_terminate = 0;
  
  // Initialize paths
  arma::mat alpha_path(q,L*L); // matrix to hold the intercepts
  arma::cube B_path(p,q,L*L);
  arma::cube Omega_path(q,q,L*L);
  arma::cube R_path(n,q,L*L);
  arma::cube tXR_path(p,q,L*L);
  arma::cube S_path(q,q,L*L);
  arma::vec theta_path(L*L);
  arma::vec eta_path(L*L);
  arma::vec obj_path(L*L); // vector holding the value of the objective function
  arma::vec prop_dir(L*L); // vector that tells us which proposal is propogated: 1 for left, 2 for down, 3 for diag
  arma::vec conv = zeros<vec>(L*L); // vector recording whether we converged or not
  arma::vec early_term = zeros<vec>(L*L); // vector recording whether we terminated early or not
  arma::vec obj_term = zeros<vec>(L*L); // vector recording whether we terminated early because objective hadn't increased substantially in consecutive iterations
  
  // First need to solve with s = 0, t = 0
  
  time_t tp;
  int time_start = time(&tp);
  s = 0;
  t = 0;
  if(verbose == 1) Rcout << "Starting s = " << s << " t = " << t << endl;
  lambda0 = lambda_spike(s);
  xi0 = xi_spike(t);
  current_index = s + L*t;
  prop_dir(current_index) = 0;
  iter = 0;
  converged = 0;
  early_terminate = 0;
  while(iter < max_iter){
    iter++;
    if(iter%50 == 0) Rcpp::checkUserInterrupt();
    B_old = B;
    Omega_old = Omega;
    old_obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty,theta_hyper_params, eta_hyper_params);
    
    // E Step
    for(int k = 0; k < q; k++){
      for(int kk = k+1; kk < q; kk++){
        tmp = 1.0/(1.0 + (1.0 - eta)/eta * xi0/xi1 * exp(-1.0 * abs(Omega(k,kk)) * (xi0 - xi1)));
        q_star(k,kk) = tmp;
        q_star(kk,k) = tmp;
      }
    }
    q_star.diag().zeros(); // to update eta, we need sum of q_star's easier to set diagonals equal to 0 and sum over entire matrix
    xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
    if(diag_penalty == 1) xi_star.diag().fill(xi1);
    else if(diag_penalty == 0) xi_star.diag().fill(0);
    
    // M Step Update of B and Theta
    update_B_theta(n,p,q,B,R,tXR,S,theta,Omega,eta,X,tXX,lambda1,lambda0,xi1,xi0,diag_penalty,theta_hyper_params,eta_hyper_params,max_iter,eps,verbose);
    // M Step Update of eta
    eta = (a_eta - 1 + accu(q_star)/2)/(a_eta + b_eta - 2 + q*(q-1)/2);
    // M Step Update of Omega
    xi_star /= n; // QUIC works with re-scaled objective and we scale our penalty accordingly
    res_quic = my_quic(q, S, xi_star, eps, quic_max_iter);
    Omega = res_quic.slice(0);
    Sigma = res_quic.slice(1);
    
    // check convergence and whether we need to terminate early
    converged = 1;
    early_terminate = 0;
    omega_non_zero = 0;
    
    obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
    s_eval = eig_sym(S);
    if(s_eval(q-1)/s_eval(0) > s_max_condition) early_terminate = 1; // condition number of S is too large. Terminate early
    if( (obj - old_obj)/old_obj < eps) obj_counter++; // objective increased by less than 100*eps%, increment counter
    else obj_counter = 0; // if objective increases by more than 100*eps% re-set counter
    
    for(int j = 0; j < p; j++){
      for(int k = 0; k < q; k++){
        if( (B_old(j,k) == 0) & (B(j,k) != 0) ) converged = 0 ; // variable was added
        else if( (B_old(j,k) != 0) & ( abs((B_old(j,k) - B(j,k))/B_old(j,k)) > eps )) converged = 0; // value of beta_{j,k} changed by more than 100*eps%
      }
    }
    for(int k = 0; k < q; k++){
      for(int kk = k+1; kk < q; kk++){
        if( (Omega_old(k,kk) == 0) & (Omega(k,kk) != 0) ) converged = 0;
        else if( (Omega_old(k,kk) != 0) & ( abs( (Omega_old(k,kk) - Omega(k,kk))/Omega_old(k,kk)) > eps)) converged = 0;
        if(Omega(k,kk) != 0) omega_non_zero++;
      }
    }
    if(obj_counter == obj_counter_max) break;
    if(early_terminate == 1) break;
    if(converged == 1) break;
  } // closes while loop of main EM algorithm
  if(obj_counter == obj_counter_max) obj_term(current_index) = 1;
  early_term(current_index) = early_terminate;
  conv(current_index) = converged;
  
  if(verbose == 1){
    if(early_terminate == 1) Rcout << "    Terminated early due to ill-conditioned S." << endl;
    if(obj_counter == obj_counter_max) Rcout << "    Terminated early due to marginal objective increase." << endl;
    if( (iter == max_iter) & (converged == 0)) Rcout <<"    Reached max iterations but failed to converge." << endl;
    Rcout << "    num B != 0 : " << accu(B != 0) << " num Omega != 0 : " << omega_non_zero;
    Rcout << "Finished s = " << s << "t = " << t << endl;
  }
  // save to paths
  B_path.slice(current_index) = B;
  Omega_path.slice(current_index) = Omega;
  R_path.slice(current_index) = R;
  tXR_path.slice(current_index) = tXR;
  S_path.slice(current_index) = S;
  theta_path(current_index) = theta;
  eta_path(current_index) = eta;
  obj_path(current_index) = objective(n, p, q, S, B, Omega, lambda1, lambda_spike(L-1), xi1, xi_spike(L-1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
  
  // Now we need to solve for t = 0
  for(s = 1; s < L; s++){
    lambda0 = lambda_spike(s);
    xi0 = xi_spike(t);
    if(verbose == 1) Rcout << "Started s = " << s << " t = " << t << endl;
    current_index = s + L*t;
    left_index = s-1 + L*t;
    if(early_term(left_index) == 1){
      // Xi^(s-1,t) is unstable. Re-start from 0
      B = B_restart;
      Omega = Omega_restart;
      R = R_restart;
      tXR = tXR_restart;
      S = S_restart;
      theta = theta_restart;
      eta = eta_restart;
    } else {
      B = B_path.slice(left_index);
      Omega = Omega_path.slice(left_index);
      R = R_path.slice(left_index);
      tXR = tXR_path.slice(left_index);
      S = S_path.slice(left_index);
      theta = theta_path(left_index);
      eta = eta_path(left_index);
    }
    prop_dir(current_index) = 1;
    iter = 0;
    converged = 0;
    early_terminate = 0;
    while(iter < max_iter){
      if(iter%50 == 0) Rcpp::checkUserInterrupt();
      iter++;
      B_old = B;
      Omega_old = Omega;
      old_obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
      // E Step
      for(int k = 0; k < q; k++){
        for(int kk = k+1; kk < q; kk++){
          tmp = 1.0/(1.0 + (1.0 - eta)/eta * xi0/xi1 * exp(-1.0 * abs(Omega(k,kk)) * (xi0 - xi1)));
          q_star(k,kk) = tmp;
          q_star(kk,k) = tmp;
        }
      }
      q_star.diag().zeros(); // to update eta, we need sum of q_star's easier to set diagonals equal to 0 and sum over entire matrix
      xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
      if(diag_penalty == 1) xi_star.diag().fill(xi1);
      else if(diag_penalty == 0) xi_star.diag().fill(0);
      
      // M Step Update of B and Theta
      update_B_theta(n,p,q,B,R,tXR,S,theta,Omega,eta,X,tXX,lambda1,lambda0,xi1,xi0,diag_penalty,theta_hyper_params,eta_hyper_params,max_iter,eps,verbose);
      // M Step Update of eta
      eta = (a_eta - 1 + accu(q_star)/2)/(a_eta + b_eta - 2 + q*(q-1)/2);
      // M Step Update of Omega
      xi_star /= n; // QUIC works with re-scaled objective and we scale our penalty accordingly
      res_quic = my_quic(q, S, xi_star, eps, quic_max_iter);
      Omega = res_quic.slice(0);
      Sigma = res_quic.slice(1);
      
      // check convergence and whether we need to terminate early
      converged = 1;
      early_terminate = 0;
      omega_non_zero = 0;
      
      obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
      s_eval = eig_sym(S);
      if(s_eval(q-1)/s_eval(0) > s_max_condition) early_terminate = 1; // condition number of S is too large. Terminate early
      if( (obj - old_obj)/old_obj < eps) obj_counter++; // objective increased by less than 100*eps%, increment counter
      else obj_counter = 0; // if objective increases by more than 100*eps% re-set counter
      
      for(int j = 0; j < p; j++){
        for(int k = 0; k < q; k++){
          if( (B_old(j,k) == 0) & (B(j,k) != 0) ) converged = 0 ; // variable was added
          else if( (B_old(j,k) != 0) & ( abs((B_old(j,k) - B(j,k))/B_old(j,k)) > eps )) converged = 0; // value of beta_{j,k} changed by more than 100*eps%
        }
      }
      for(int k = 0; k < q; k++){
        for(int kk = k+1; kk < q; kk++){
          if( (Omega_old(k,kk) == 0) & (Omega(k,kk) != 0) ) converged = 0;
          else if( (Omega_old(k,kk) != 0) & ( abs( (Omega_old(k,kk) - Omega(k,kk))/Omega_old(k,kk)) > eps)) converged = 0;
          if(Omega(k,kk) != 0) omega_non_zero++;
        }
      }
      if(obj_counter == obj_counter_max) break;
      if(early_terminate == 1) break;
      if(converged == 1) break;
    } // closes while loop of main EM algorithm
    if(obj_counter == obj_counter_max) obj_term(current_index) = 1;
    early_term(current_index) = early_terminate;
    conv(current_index) = converged;
    if(verbose == 1){
      if(early_terminate == 1) Rcout << "    Terminated early due to ill-conditioned S." << endl;
      if(obj_counter == obj_counter_max) Rcout << "    Terminated early due to marginal objective increase." << endl;
      if((iter == max_iter) & (converged == 0)) Rcout << "    Reached max iterations but failed to converge." << endl;
      Rcout << "    num B != 0 : " << accu(B != 0) << " num Omega != 0 : " << omega_non_zero;
      Rcout << "Finished s = " << s << "t = " << t << endl;
    }
    // save to paths
    B_path.slice(current_index) = B;
    Omega_path.slice(current_index) = Omega;
    R_path.slice(current_index) = R;
    tXR_path.slice(current_index) = tXR;
    S_path.slice(current_index) = S;
    theta_path(current_index) = theta;
    eta_path(current_index) = eta;
    obj_path(current_index) = objective(n, p, q, S, B, Omega, lambda1, lambda_spike(L-1), xi1, xi_spike(L-1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
    
  }
  // Now we need to solve for s = 0
  s = 0;
  t = 0;
  for(t = 1; t < L; t++){
    lambda0 = lambda_spike(s);
    xi0 = xi_spike(t);
    if(verbose == 1) Rcout << "Started s = " << s << " t = " << t << endl;
    current_index = s + L*t;
    down_index = s + L*(t-1);
    if(early_term(down_index) == 1){
      B = B_restart;
      Omega = Omega_restart;
      R = R_restart;
      tXR = tXR_restart;
      S = S_restart;
      theta = theta_restart;
      eta = eta_restart;
    } else {
      B = B_path.slice(down_index);
      Omega = Omega_path.slice(down_index);
      R = R_path.slice(down_index);
      tXR = tXR_path.slice(down_index);
      S = S_path.slice(down_index);
      theta = theta_path(down_index);
      eta = eta_path(down_index);
    }
    prop_dir(current_index) = 2;
    iter = 0;
    converged = 0;
    early_terminate = 0;
    while(iter < max_iter){
      iter++;
      if(iter%50) Rcpp::checkUserInterrupt();
      B_old = B;
      Omega_old = Omega;
      old_obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
      
      // E Step
      for(int k = 0; k < q; k++){
        for(int kk = k+1; kk < q; kk++){
          tmp = 1.0/(1.0 + (1.0 - eta)/eta * xi0/xi1 * exp(-1.0 * abs(Omega(k,kk)) * (xi0 - xi1)));
          q_star(k,kk) = tmp;
          q_star(kk,k) = tmp;
        }
      }
      q_star.diag().zeros(); // to update eta, we need sum of q_star's easier to set diagonals equal to 0 and sum over entire matrix
      xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
      if(diag_penalty == 1) xi_star.diag().fill(xi1);
      else if(diag_penalty == 0) xi_star.diag().fill(0);
      
      // M Step Update of B and Theta
      update_B_theta(n,p,q,B,R,tXR,S,theta,Omega,eta,X,tXX,lambda1,lambda0,xi1,xi0,diag_penalty,theta_hyper_params,eta_hyper_params,max_iter,eps,verbose);
      // M Step Update of eta
      eta = (a_eta - 1 + accu(q_star)/2)/(a_eta + b_eta - 2 + q*(q-1)/2);
      // M Step Update of Omega
      xi_star /= n; // QUIC works with re-scaled objective and we scale our penalty accordingly
      res_quic = my_quic(q, S, xi_star, eps, quic_max_iter);
      Omega = res_quic.slice(0);
      Sigma = res_quic.slice(1);
      
      // check convergence and whether we need to terminate early
      converged = 1;
      early_terminate = 0;
      omega_non_zero = 0;
      
      obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
      s_eval = eig_sym(S);
      if(s_eval(q-1)/s_eval(0) > s_max_condition) early_terminate = 1; // condition number of S is too large. Terminate early
      if( (obj - old_obj)/old_obj < eps) obj_counter++; // objective increased by less than 100*eps%, increment counter
      else obj_counter = 0; // if objective increases by more than 100*eps% re-set counter
      
      for(int j = 0; j < p; j++){
        for(int k = 0; k < q; k++){
          if( (B_old(j,k) == 0) & (B(j,k) != 0) ) converged = 0 ; // variable was added
          else if( (B_old(j,k) != 0) & ( abs((B_old(j,k) - B(j,k))/B_old(j,k)) > eps )) converged = 0; // value of beta_{j,k} changed by more than 100*eps%
        }
      }
      for(int k = 0; k < q; k++){
        for(int kk = k+1; kk < q; kk++){
          if( (Omega_old(k,kk) == 0) & (Omega(k,kk) != 0) ) converged = 0;
          else if( (Omega_old(k,kk) != 0) & ( abs( (Omega_old(k,kk) - Omega(k,kk))/Omega_old(k,kk)) > eps)) converged = 0;
          if(Omega(k,kk) != 0) omega_non_zero++;
        }
      }
      if(obj_counter == obj_counter_max) break;
      if(early_terminate == 1) break;
      if(converged == 1) break;
    } // closes while loop of main EM algorithm
    if(obj_counter == obj_counter_max) obj_term(current_index) = 1;
    early_term(current_index) = early_terminate;
    conv(current_index) = converged;
    if(verbose == 1){
      if(early_terminate == 1) Rcout << "    Terminated early due to ill-conditioned S." << endl;
      if(obj_counter == obj_counter_max) Rcout << "    Terminated early due to marginal objective increase." << endl;
      if((iter == max_iter) & (converged == 0)) Rcout << "    Reached max iterations but failed to converge." << endl;
      Rcout << "    num B != 0 : " << accu(B != 0) << " num Omega != 0 : " << omega_non_zero;
      Rcout << "Finished s = " << s << "t = " << t << endl;
    }
    // save to paths
    B_path.slice(current_index) = B;
    Omega_path.slice(current_index) = Omega;
    R_path.slice(current_index) = R;
    tXR_path.slice(current_index) = tXR;
    S_path.slice(current_index) = S;
    theta_path(current_index) = theta;
    eta_path(current_index) = eta;
    obj_path(current_index) = objective(n, p, q, S, B, Omega, lambda1, lambda_spike(L-1), xi1, xi_spike(L-1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
    
  }
  // Now we're ready to solve with s > 0 and t > 0
  
  for(s = 1; s < L; s++){
    for(t = 1; t < L; t++){
      if(verbose == 1) Rcout << "Started s = " << s << " t = " << t << endl;
      lambda0 = lambda_spike(s);
      xi0 = xi_spike(t);
      current_index = s + L*t;
      
      left_index = s-1 + L*t;
      down_index = s + L*(t-1);
      diag_index = (s-1) + L*(t-1);
      
      if(early_term(left_index) == 1){ // solution to left was unstable. use the re-start
        B_left = B_restart;
        Omega_left = Omega_restart;
        R_left = R_restart;
        tXR_left = tXR_restart;
        S_left = S_restart;
        theta_left = theta_restart;
        eta_left = eta_restart;
      } else{
        B_left = B_path.slice(left_index);
        Omega_left = Omega_path.slice(left_index);
        R_left = R_path.slice(left_index);
        tXR_left = tXR_path.slice(left_index);
        S_left = S_path.slice(left_index);
        theta_left = theta_path(left_index);
        eta_left = eta_path(left_index);
      }
      obj_left = objective(n, p, q, S_left, B_left, Omega_left, lambda1, lambda0, xi1, xi0, theta_left, eta_left, diag_penalty, theta_hyper_params, eta_hyper_params);
      if(early_term(down_index) == 1){
        B_down = B_restart;
        Omega_down = Omega_restart;
        R_down = R_restart;
        tXR_down = tXR_restart;
        S_down = S_restart;
        theta_down = theta_restart;
        eta_down = eta_restart;
      } else {
        B_down = B_path.slice(down_index);
        Omega_down = Omega_path.slice(down_index);
        R_down = R_path.slice(down_index);
        tXR_down = tXR_path.slice(down_index);
        S_down = S_path.slice(down_index);
        theta_down = theta_path(down_index);
        eta_down = eta_path(down_index);
      }
      obj_down = objective(n,p,q,S_down, B_down, Omega_down, lambda1, lambda0, xi1, xi0, theta_down, eta_down, diag_penalty, theta_hyper_params, eta_hyper_params);
      if(early_term(diag_index) == 1){
        B_diag = B_restart;
        Omega_diag = Omega_restart;
        R_diag = R_restart;
        tXR_diag = tXR_restart;
        S_diag = S_restart;
        theta_diag = theta_restart;
        eta_diag = eta_restart;
      } else {
        B_diag = B_path.slice(diag_index);
        Omega_diag = Omega_path.slice(diag_index);
        R_diag = R_path.slice(diag_index);
        tXR_diag = tXR_path.slice(diag_index);
        S_diag = S_path.slice(diag_index);
        theta_diag = theta_path(diag_index);
        eta_diag = eta_path(diag_index);
      }
      obj_diag = objective(n, p, q, S_diag, B_diag, Omega_diag, lambda1, lambda0, xi1, xi0, theta_diag, eta_diag, diag_penalty, theta_hyper_params, eta_hyper_params);
      
      if( (obj_left >= obj_down) & (obj_left >= obj_diag) ){
        B = B_left;
        Omega = Omega_left;
        R = R_left;
        tXR = tXR_left;
        S = S_left;
        theta = theta_left;
        eta = eta_left;
        prop_dir(current_index) = 1;
        old_obj = obj_left;
      } else if( (obj_down > obj_left) & (obj_down >= obj_diag) ){
        B = B_down;
        Omega = Omega_down;
        R = R_down;
        tXR = tXR_down;
        S = S_down;
        theta = theta_down;
        eta = eta_down;
        prop_dir(current_index) = 2;
        old_obj = obj_down;
      } else if( (obj_diag > obj_left) & (obj_diag > obj_down) ){
        B = B_diag;
        Omega = Omega_diag;
        R = R_diag;
        tXR = tXR_diag;
        S = S_diag;
        theta = theta_diag;
        eta = eta_diag;
        prop_dir(current_index) = 3;
        old_obj = obj_diag;
      }
      iter = 0;
      converged = 0;
      early_terminate = 0;
      while(iter < max_iter){
        iter++;
        if(iter%50) Rcpp::checkUserInterrupt();
        B_old = B;
        Omega_old = Omega;
        old_obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
        
        // E Step
        for(int k = 0; k < q; k++){
          for(int kk = k+1; kk < q; kk++){
            tmp = 1.0/(1.0 + (1.0 - eta)/eta * xi0/xi1 * exp(-1.0 * abs(Omega(k,kk)) * (xi0 - xi1)));
            q_star(k,kk) = tmp;
            q_star(kk,k) = tmp;
          }
        }
        q_star.diag().zeros(); // to update eta, we need sum of q_star's easier to set diagonals equal to 0 and sum over entire matrix
        xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
        if(diag_penalty == 1) xi_star.diag().fill(xi1);
        else if(diag_penalty == 0) xi_star.diag().fill(0);
        
        // M Step Update of B and Theta
        update_B_theta(n,p,q,B,R,tXR,S,theta,Omega,eta,X,tXX,lambda1,lambda0,xi1,xi0,diag_penalty,theta_hyper_params,eta_hyper_params,max_iter,eps,verbose);
        // M Step Update of eta
        eta = (a_eta - 1 + accu(q_star)/2)/(a_eta + b_eta - 2 + q*(q-1)/2);
        // M Step Update of Omega
        xi_star /= n; // QUIC works with re-scaled objective and we scale our penalty accordingly
        res_quic = my_quic(q, S, xi_star, eps, quic_max_iter);
        Omega = res_quic.slice(0);
        Sigma = res_quic.slice(1);
        
        // check convergence and whether we need to terminate early
        converged = 1;
        early_terminate = 0;
        omega_non_zero = 0;
        
        obj = objective(n, p, q, S, B, Omega, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
        s_eval = eig_sym(S);
        if(s_eval(q-1)/s_eval(0) > s_max_condition) early_terminate = 1; // condition number of S is too large. Terminate early
        if( (obj - old_obj)/old_obj < eps) obj_counter++; // objective increased by less than 100*eps%, increment counter
        else obj_counter = 0; // if objective increases by more than 100*eps% re-set counter
        
        for(int j = 0; j < p; j++){
          for(int k = 0; k < q; k++){
            if( (B_old(j,k) == 0) & (B(j,k) != 0) ) converged = 0 ; // variable was added
            else if( (B_old(j,k) != 0) & ( abs((B_old(j,k) - B(j,k))/B_old(j,k)) > eps )) converged = 0; // value of beta_{j,k} changed by more than 100*eps%
          }
        }
        for(int k = 0; k < q; k++){
          for(int kk = k+1; kk < q; kk++){
            if( (Omega_old(k,kk) == 0) & (Omega(k,kk) != 0) ) converged = 0;
            else if( (Omega_old(k,kk) != 0) & ( abs( (Omega_old(k,kk) - Omega(k,kk))/Omega_old(k,kk)) > eps)) converged = 0;
            if(Omega(k,kk) != 0) omega_non_zero++;
          }
        }
        if(obj_counter == obj_counter_max) break;
        if(early_terminate == 1) break;
        if(converged == 1) break;
      } // closes while loop of main EM algorithm
      if(obj_counter == obj_counter_max) obj_term(current_index) = 1;
      early_term(current_index) = early_terminate;
      conv(current_index) = converged;
      if(verbose == 1){
        if(early_terminate == 1) Rcout << "    Terminated early due to ill-conditioned S." << endl;
        if(obj_counter == obj_counter_max) Rcout << "    Terminated early due to marginal objective increase." << endl;
        if((iter == max_iter) & (converged == 0)) Rcout << "    Reached max iterations but failed to converge." << endl;
        Rcout << "    num B != 0 : " << accu(B != 0) << " num Omega != 0 : " << omega_non_zero;
        Rcout << "Finished s = " << s << "t = " << t << endl;
      }
      // save to paths
      B_path.slice(current_index) = B;
      Omega_path.slice(current_index) = Omega;
      R_path.slice(current_index) = R;
      tXR_path.slice(current_index) = tXR;
      S_path.slice(current_index) = S;
      theta_path(current_index) = theta;
      eta_path(current_index) = eta;
      obj_path(current_index) = objective(n, p, q, S, B, Omega, lambda1, lambda_spike(L-1), xi1, xi_spike(L-1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params); // compute at the maximal value
    } // closes the loop over t
  } // closes the loop over s
  
  int time_end = time(&tp);
  List results;
  for(int s = 0; s < L; s++){
    for(int t = 0; t < L; t++){
      current_index = s + L*t;
      tmp_B = B_path.slice(current_index);
      tmp_B.each_col()/= x_col_weights;
      B_path.slice(current_index) = tmp_B;
      alpha_path.col(current_index) = mu_y - tmp_B.t()*mu_x; // compute the intercept
    }
  }
  alpha = mu_y - tmp_B.t() * mu_x;
  results["alpha"] = alpha;
  results["B"] = B_path.slice(L*L - 1);
  results["Omega"] = Omega_path.slice(L*L-1);
  //results["R"] = R_path.slice(L*L-1);
  //results["S"] = S_path.slice(L*L-1);
  results["theta"] = theta_path(L*L-1);
  results["eta"] = eta_path(L*L-1);
  results["alpha_path"] = alpha_path;
  results["B_path"] = B_path;
  results["Omega_path"] = Omega_path;
  //results["R_path"] = R_path;
  //results["S_path"] = S_path;
  results["theta_path"] = theta_path;
  results["eta_path"] = eta_path;
  //results["obj_path"] = obj_path;
  //results["prop_dir"] = prop_dir;
  //results["conv"] = conv;
  results["time"] = time_end - time_start;
  results["lambda0"] = lambda_spike;
  results["xi0"] = xi_spike;
  results["early_term"] = early_term;
  //results["obj_term"] = obj_term;
  return(results);
}

// [[Rcpp::export]]
List mSSL_dcpe(arma::mat X,
               arma::mat Y,
               List lambdas,
               List xis,
               arma::vec theta_hyper_params,
               arma::vec eta_hyper_params,
               int diag_penalty,
               int max_iter,
               double eps,
               int verbose)
{

  int n = Y.n_rows;
  int q = Y.n_cols;
  int p = X.n_cols;
  
  // Center columns of X and Y
  // re-scale the re-centered columns of X to have norm sqrt(n)
  
  double tmp_mu_x = 0.0;
  double tmp_weight_x = 0.0;
  double tmp_mu_y = 0.0;
  arma::vec x_col_weights(p);
  arma::vec mu_x(p);
  arma::vec mu_y(q);
  for(int j = 0; j < p ; j++){
    tmp_mu_x = mean(X.col(j));
    X.col(j) -= tmp_mu_x;
    tmp_weight_x = norm(X.col(j))/sqrt(n);
    X.col(j) /= tmp_weight_x;
    mu_x(j) = tmp_mu_x;
    x_col_weights(j) = tmp_weight_x;
  }
  for(int k = 0; k < q; k++){
    tmp_mu_y = mean(Y.col(k));
    Y.col(k) -= tmp_mu_y;
    mu_y(k) = tmp_mu_y;
  }
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
  
  arma::mat tYY = Y.t() * Y;
  arma::mat tXX = X.t() * X;
  arma::mat tXY = X.t() * Y;
  arma::mat R = Y - X * B;
  arma::mat tXR = X.t() * R;
  arma::mat tRR = R.t() * R;
  arma::mat S = tRR/n;
  if(verbose == 1) Rcout << "Initialized R, tXR, S" << endl;
  
  
  // if we ever fail to converge, we should re-set the values of B and Omega
  arma::mat B_reset = B;
  arma::mat Omega_reset = Omega;
  arma::mat Sigma_reset = Sigma;
  arma::mat R_reset = R;
  arma::mat tXR_reset = tXR;
  arma::mat tRR_reset = tRR;
  arma::mat S_reset = S;
  double theta_reset = theta;
  double eta_reset = eta;
  
  //Rcout << "preparing QUIC working parameters" << endl;
  // initialize stuff for QUIC
  cube res_quic(q,q,2); // cube to hold values from QUIC. 1st one is for Omega, 2nd is for Sigma
  //double* S_ptr = S.memptr();
  //double* Xi_ptr = xi_star.memptr();
  
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
    if(verbose == 1) Rcout << setprecision(6) << "Starting lambda0 = " << lambda0 << " num B non-zero = " << accu(B != 0) << " theta = " << theta << endl;
    update_B_theta(n,p,q,B,R,tXR,S,theta,Omega,eta,X,tXX,lambda1,lambda0,xi1,xi0,diag_penalty,theta_hyper_params,eta_hyper_params,max_iter,eps, verbose);
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
      eta = (a_eta - 1 + accu(q_star)/2)/(a_eta + b_eta -2 + q*(q-1)/2);
      // M-step update of Omega
      xi_star /= n; // QUIC needs everything to be scaled
      res_quic = my_quic(q, S, xi_star, eps, quic_max_iter);
      Omega = res_quic.slice(0);
      Sigma = res_quic.slice(1);
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
  R_reset = R;
  tXR_reset = tXR;
  S_reset = S;
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
    
    update_B_theta(n,p,q,B,R,tXR,S,theta,Omega,eta,X,tXX,lambda1,lambda0,xi1,xi0,diag_penalty,theta_hyper_params,eta_hyper_params,max_iter,eps, verbose);
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
    eta = (a_eta - 1 + accu(q_star)/2)/(a_eta + b_eta -2 + q*(q-1)/2);
    // M-step update of Omega
    xi_star /= n;
    res_quic = my_quic(q, S, xi_star, eps, quic_max_iter);
    Omega = res_quic.slice(0);
    Sigma = res_quic.slice(1);
    
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
    if(verbose == 1) Rcout << "iter = " << iter << " num nonzero B: " << accu(B!=0) << " num nonzero Omega: " << omega_non_zero << " theta = " << theta << " eta = " << eta << endl;
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
  alpha = mu_y - tmp_B.t() * mu_x;
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


// [[Rcpp::export]]
List gSSL(arma::mat Y,
          List xis,
          arma::vec eta_hyper_params,
          int diag_penalty,
          int max_iter,
          double eps,
          int obj_counter_max,
          int verbose){
  
  int n = Y.n_rows;
  int q = Y.n_cols;
  int quic_max_iter = 5*max_iter;

  /*
  // Read in the control parameters
  int max_iter = control_params["max_iter"];
  double eps = control_params["eps"];
  int obj_counter_max = control_params["obj_counter_max"];
  */
  
  // Initialize grid of xi0 values
  double xi1 = xis["xi1"];
  arma::vec xi_spike = xis["xi0"];
  double xi0 = xi_spike(0);
  int L = xi_spike.n_elem;
  
  // Read hyper-parameters
  double a_eta = eta_hyper_params(0);
  double b_eta = eta_hyper_params(1);
  
  // Initialize some working parameters
  arma::mat Omega(q,q);
  Omega.eye();
  arma::mat Sigma(q,q);
  Sigma.eye();
  arma::mat Omega_old(q,q);
  Omega_old.eye();
  
  double eta = 0.5;
  double obj = 0.0; // value of objective function
  double obj_old = 0.0; // used to check convergence
  
  arma::mat tYY = Y.t() * Y;
  arma::mat S = tYY/n;
  arma::vec s_eval = eig_sym(S); // vector of eigenvalues of S
  
  // Initialize q_star and penalty matrices
  arma::mat q_star(q,q);
  q_star.fill(0.5);
  arma::mat xi_star(q,q);
  xi_star.fill(xi1/n);
  xi_star.diag().zeros(); // we will eventually penalize the diagonal
  double tmp = 0.0; // holds the value of q_star(k,k');
  int omega_non_zero = 0; // counter for number of non-zero off-diagonal elements in Omega
  
  // Initialize stuff for QUIC
  cube res_quic(q,q,2);
  
  
  // Initialize iterators and counters
  int t = 0;
  int iter = 0;
  int converged = 0;
  int obj_counter = 0;
  
  // Initialize paths
  arma::cube Omega_path(q,q,L);
  arma::cube Sigma_path(q,q,L);
  arma::vec eta_path(L);
  
  time_t tp;
  int time_start = time(&tp);
  for(t = 0; t < L; t++){
    xi0 = xi_spike(t);
    if(verbose == 1) Rcout << "Starting xi0 = " << xi0 << endl;
    iter = 0;
    converged = 0;
    obj_counter = 0;
    while(iter < max_iter){
      iter++;
      Omega_old = Omega;
      obj_old = gSSL_objective(n, q, S, Omega, xi1, xi0, eta, diag_penalty, eta_hyper_params);
      
      // E Step
      for(int k = 0; k < q; k++){
        for(int kk = k+1; kk < q; kk++){
          tmp = 1.0/(1.0 + (1.0 - eta)/eta * xi0/xi1 * exp(-1.0 * abs(Omega(k,kk)) * (xi0  - xi1)));
          q_star(k,kk) = tmp;
          q_star(kk,k) = tmp;
        }
      }
      q_star.diag().zeros(); // to update eta, we need the sum of q_star's easier to set diagonals equal to 0 and usm over entire matrix
      xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
      if(diag_penalty == 1) xi_star.diag().fill(xi1);
      else if(diag_penalty == 0) xi_star.diag().fill(0);
      
      // M Step
      eta = (a_eta - 1 + accu(q_star)/2)/(a_eta + b_eta - 2 + q*(q-1)/2);
      xi_star /= n;
      res_quic = my_quic(q, S, xi_star, eps, quic_max_iter);
      Omega = res_quic.slice(0);
      Sigma = res_quic.slice(1);
      
      // check convergence
      converged = 1;
      omega_non_zero = 0;
      
      obj = gSSL_objective(n, q, S, Omega, xi1, xi0, eta, diag_penalty, eta_hyper_params);
      if( (obj - obj_old)/obj_old < eps ) obj_counter++; //objective increased by less than 100*eps, incremement counter.
      else obj_counter = 0;
      
      for(int k = 0; k < q; k++){
        for(int kk = k+1; kk < q; kk++){
          if( (Omega_old(k,kk) == 0) & (Omega(k,kk) != 0) ) converged = 0;
          else if( (Omega_old(k,kk) != 0) & ( abs((Omega_old(k,kk) - Omega(k,kk))/Omega_old(k,kk)) > eps)) converged = 0;
          if(Omega(k,kk) != 0) omega_non_zero++;
        }
      }
      if(obj_counter == obj_counter_max) break;
      if(converged == 1) break;
    } // Closes main EM loop
    if(verbose == 1){
      if(obj_counter == obj_counter_max) Rcout << "    Terminated early due to marginal objective increase." << endl;
      if( (iter == max_iter) & (converged == 0) ) Rcout << "    Reached max iteration but failed to converged" << endl;
      Rcout << "  num Omega != 0 : " << omega_non_zero << endl;
    }
    // save to paths
    Omega_path.slice(t) = Omega;
    Sigma_path.slice(t) = Sigma;
    eta_path(t) = eta;
  } // close loop over xi0's
  int time_end = time(&tp);
  List results;
  results["Omega"] = Omega;
  results["Sigma"] = Sigma;
  results["eta"] = eta;
  results["Omega_path"] = Omega_path;
  results["Sigma_path"] = Sigma_path;
  results["eta_path"] = eta_path;
  results["time"] = time_end - time_start;
  return(results);
}
