
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
List gSSLcpp(arma::mat Y,
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
