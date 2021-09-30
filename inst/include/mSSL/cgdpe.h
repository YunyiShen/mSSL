#ifndef MSSL_CGDPE_H
#define MSSL_CGDPE_H

#include <RcppArmadillo.h>
#include <linconGaussR.h>
#include <QUIC.h>
#include <workingparam.h>
#include "B_theta_update.h"

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace workingparam;
using namespace linconGaussR;
using namespace workingparam;

namespace mSSL{

template<class workpara>
List cgSSL_dpe(workpara instance,
               List lambdas,
               List xis,
               arma::vec theta_hyper_params,
               arma::vec eta_hyper_params,
               int diag_penalty,
               int max_iter,
               double eps,
               int s_max_condition,
               int obj_counter_max,
               int verbose,
               int n_rep = 1000,
               int nskp = 5)
{
  int n = instance.n_Omega;
  int q = instance.Y.n_cols;
  int p = instance.X.n_cols;
  
  // Center columns of X and Y
  // re-scale the re-centered columns of X to have norm sqrt(n)

  double tmp_mu_x = 0.0;
  double tmp_weight_x = 0.0;
  //double tmp_mu_y = 0.0;
  arma::vec x_col_weights(p);
  arma::vec mu_x(p);
  for (int j = 0; j < p; j++)
  {
    tmp_mu_x = mean(instance.X.col(j));
    instance.X.col(j) -= tmp_mu_x;
    tmp_weight_x = norm(instance.X.col(j)) / sqrt(instance.n_B);
    instance.X.col(j) /= tmp_weight_x;
    mu_x(j) = tmp_mu_x;
    x_col_weights(j) = tmp_weight_x;
  }
  instance.tXX = instance.X.t() * instance.X;
  int quic_max_iter = 5 * max_iter; // allow QUIC to run for more iterations

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
  arma::mat B(p, q);
  B.zeros();
  arma::mat B_old = B; // used internally to check convergence
  arma::mat tmp_B = B; // used at the end to re-scale columns
  arma::mat Omega(q, q);
  Omega.eye();
  arma::mat Omega_old = Omega;
  arma::mat Sigma(q, q);
  Sigma.eye();
  double theta = 0.5;
  double eta = 0.5;
  double obj = 0.0;     // value of objective function
  double old_obj = 0.0; // used to check convergence
  arma::vec mu_old(q);

  
  // Initialize a copy for the re-starts
  arma::mat B_restart = B;
  arma::mat Omega_restart = Omega;
  arma::mat Sigma_restart = Sigma;
  workpara instance_restart = instance;
  double theta_restart = theta;
  double eta_restart = eta;

  // Initialize p_star and q_star and penalty matrices
  arma::mat q_star(q, q);
  q_star.fill(0.5);
  arma::mat xi_star(q, q);
  xi_star.fill(xi1 / n);
  xi_star.diag().zeros();
  double tmp = 0.0;       // holds the value of q_star(k,k')
  int omega_non_zero = 0; // counter for number of non-zero off-diagonal entries in Omega

  // Initialize stuff for QUIC
  cube res_quic(q, q, 2); // cube to hold value of Omega and Sigma computed in QUIC

  // We need copies of B, Omega, R, tXR, S, theta, and eta that are propogated along the grid
  arma::mat B_left = B;
  arma::mat Omega_left = Omega;
  arma::mat Sigma_left = Sigma;
  workpara instance_left = instance;
  double theta_left = theta;
  double eta_left = eta;
  int left_index = 0;
  double obj_left = 0.0;

  arma::mat B_down = B;
  arma::mat Omega_down = Omega;
  arma::mat Sigma_down = Sigma;
  workpara instance_down = instance;
  double theta_down = theta;
  double eta_down = eta;
  int down_index = 0;
  double obj_down = 0.0;

  arma::mat B_diag = B;
  arma::mat Omega_diag = Omega;
  arma::mat Sigma_diag = Sigma;
  workpara instance_diag = instance;
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
  //  YS: not quite sure if this is needed, they have huge memory footprints
  arma::mat alpha_path(q, L * L); // matrix to hold the intercepts
  arma::cube B_path(p, q, L * L);
  arma::cube Omega_path(q, q, L * L);
  arma::cube Sigma_path(q, q, L * L);
  arma::cube R_path(instance.R.n_rows, q, L * L);
  arma::cube tXR_path(p, q, L * L);
  arma::cube S_Omega_path(q, q, L * L);
  arma::cube S_B_path(q, q, L * L);
  arma::cube M_path(q, q, L * L);
  arma::vec theta_path(L * L);
  arma::vec eta_path(L * L);
  arma::vec obj_path(L * L);                // vector holding the value of the objective function
  arma::vec prop_dir(L * L);                // vector that tells us which proposal is propogated: 1 for left, 2 for down, 3 for diag
  arma::vec conv = zeros<vec>(L * L);       // vector recording whether we converged or not
  arma::vec early_term = zeros<vec>(L * L); // vector recording whether we terminated early or not
  arma::vec obj_term = zeros<vec>(L * L);   // vector recording whether we terminated early because objective hadn't increased substantially in consecutive iterations

  // First need to solve with s = 0, t = 0

  time_t tp;
  int time_start = time(&tp);
  s = 0;
  t = 0;
  if (verbose == 1)
    Rcout << "Starting s = " << s << " t = " << t << endl;
  lambda0 = lambda_spike(s);
  xi0 = xi_spike(t);
  current_index = s + L * t;
  prop_dir(current_index) = 0;
  iter = 0;
  converged = 0;
  early_terminate = 0;
  while (iter < max_iter)
  {
    iter++;
    if (iter % 50 == 0)
      Rcpp::checkUserInterrupt();
    B_old = B;
    Omega_old = Omega;
    old_obj = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
    

    // E Step
    for (int k = 0; k < q; k++)
    {
      for (int kk = k + 1; kk < q; kk++)
      {
        tmp = 1.0 / (1.0 + (1.0 - eta) / eta * xi0 / xi1 * exp(-1.0 * abs(Omega(k, kk)) * (xi0 - xi1)));
        q_star(k, kk) = tmp;
        q_star(kk, k) = tmp;
      }
    }
    q_star.diag().zeros(); // to update eta, we need sum of q_star's easier to set diagonals equal to 0 and sum over entire matrix
    xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
    if (diag_penalty == 1)
      xi_star.diag().fill(xi1);
    else if (diag_penalty == 0)
      xi_star.diag().fill(0);

    // M Step Update of B and Theta
    update_B_theta(instance.n_B, p, q, B, instance.R, instance.tXR, instance.S, theta, Sigma, eta, instance.X, instance.tXX, lambda1, lambda0, xi1, xi0, diag_penalty, theta_hyper_params, eta_hyper_params, max_iter, eps, verbose);
    instance.update_M(B);
    // M Step Update of eta
    eta = (a_eta - 1 + arma::accu(q_star) / 2) / (a_eta + b_eta - 2 + q * (q - 1) / 2);
    // M Step Update of Omega
    xi_star /= n; // QUIC works with re-scaled objective and we scale our penalty accordingly
    res_quic = cgquic(q, instance.S_Omega, instance.M, xi_star, eps, quic_max_iter);
    Omega = res_quic.slice(0);
    Sigma = res_quic.slice(1);
    instance.postprocessing(B,Sigma,Omega);
    
    if(verbose == 1) Rcout << "    updating residual matrices "<< endl;
    instance.update(mu_old,B,Sigma,Omega,n_rep,nskp);
    mu_old = instance.mu;

    // check convergence and whether we need to terminate early
    converged = 1;
    early_terminate = 0;
    omega_non_zero = 0;

    obj = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
    
    if (instance.s_eval(q - 1) / instance.s_eval(0) > s_max_condition)
      early_terminate = 1; // condition number of S is too large. Terminate early
    if ((obj - old_obj) / old_obj < eps)
      obj_counter++; // objective increased by less than 100*eps%, increment counter
    else
      obj_counter = 0; // if objective increases by more than 100*eps% re-set counter

    for (int j = 0; j < p; j++)
    {
      for (int k = 0; k < q; k++)
      {
        if ((B_old(j, k) == 0) & (B(j, k) != 0))
          converged = 0; // variable was added
        else if ((B_old(j, k) != 0) & (abs((B_old(j, k) - B(j, k)) / B_old(j, k)) > eps))
          converged = 0; // value of beta_{j,k} changed by more than 100*eps%
      }
    }
    for (int k = 0; k < q; k++)
    {
      for (int kk = k + 1; kk < q; kk++)
      {
        if ((Omega_old(k, kk) == 0) & (Omega(k, kk) != 0))
          converged = 0;
        else if ((Omega_old(k, kk) != 0) & (abs((Omega_old(k, kk) - Omega(k, kk)) / Omega_old(k, kk)) > eps))
          converged = 0;
        if (Omega(k, kk) != 0)
          omega_non_zero++;
      }
    }
    if (obj_counter == obj_counter_max)
      break;
    if (early_terminate == 1)
      break;
    if (converged == 1)
      break;
  } // closes while loop of main EM algorithm
  if (obj_counter == obj_counter_max)
    obj_term(current_index) = 1;
  early_term(current_index) = early_terminate;
  conv(current_index) = converged;

  if (verbose == 1)
  {
    if (early_terminate == 1)
      Rcout << "    Terminated early due to ill-conditioned S." << endl;
    if (obj_counter == obj_counter_max)
      Rcout << "    Terminated early due to marginal objective increase." << endl;
    if ((iter == max_iter) & (converged == 0))
      Rcout << "    Reached max iterations but failed to converge." << endl;
    Rcout << "    num B != 0 : " << arma::accu(B != 0) << " num Omega != 0 : " << omega_non_zero;
    Rcout << "Finished s = " << s << "t = " << t << endl;
  }
  // save to paths
  B_path.slice(current_index) = B;
  Omega_path.slice(current_index) = Omega;
  Sigma_path.slice(current_index) = Sigma;
  R_path.slice(current_index) = instance.R;
  tXR_path.slice(current_index) = instance.tXR;
  S_B_path.slice(current_index) = instance.S;
  S_Omega_path.slice(current_index) = instance.S_Omega;
  M_path.slice(current_index) = instance.M;
  theta_path(current_index) = theta;
  eta_path(current_index) = eta;

  obj_path(current_index) = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda_spike(L - 1), xi1, xi_spike(L - 1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);

  // Now we need to solve for t = 0
  for (s = 1; s < L; s++)
  {
    lambda0 = lambda_spike(s);
    xi0 = xi_spike(t);
    if (verbose == 1)
      Rcout << "Started s = " << s << " t = " << t << endl;
    current_index = s + L * t;
    left_index = s - 1 + L * t;
    if (early_term(left_index) == 1)
    {
      // Xi^(s-1,t) is unstable. Re-start from 0
      B = B_restart;
      Omega = Omega_restart;
      Sigma = Omega_restart;
      instance = instance_restart;
      theta = theta_restart;
      eta = eta_restart;
    }
    else
    {
      B = B_path.slice(left_index);
      Omega = Omega_path.slice(left_index);
      Sigma = Sigma_path.slice(left_index);
      instance.R = R_path.slice(left_index);
      instance.tXR = tXR_path.slice(left_index);
      instance.S = S_B_path.slice(left_index);
      instance.S_Omega = S_Omega_path.slice(left_index);
      instance.M = M_path.slice(left_index);
      theta = theta_path(left_index);
      eta = eta_path(left_index);
    }
    prop_dir(current_index) = 1;
    iter = 0;
    converged = 0;
    early_terminate = 0;
    while (iter < max_iter)
    {
      if (iter % 50 == 0)
        Rcpp::checkUserInterrupt();
      iter++;
      B_old = B;
      Omega_old = Omega;

      old_obj = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
      // E Step
      for (int k = 0; k < q; k++)
      {
        for (int kk = k + 1; kk < q; kk++)
        {
          tmp = 1.0 / (1.0 + (1.0 - eta) / eta * xi0 / xi1 * exp(-1.0 * abs(Omega(k, kk)) * (xi0 - xi1)));
          q_star(k, kk) = tmp;
          q_star(kk, k) = tmp;
        }
      }
      q_star.diag().zeros(); // to update eta, we need sum of q_star's easier to set diagonals equal to 0 and sum over entire matrix
      xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
      if (diag_penalty == 1)
        xi_star.diag().fill(xi1);
      else if (diag_penalty == 0)
        xi_star.diag().fill(0);

      
      // M Step Update of B and Theta
      update_B_theta(instance.n_B, p, q, B, instance.R, instance.tXR, instance.S, theta, Sigma, eta, instance.X, instance.tXX, lambda1, lambda0, xi1, xi0, diag_penalty, theta_hyper_params, eta_hyper_params, max_iter, eps, verbose);
      // update things related with B
      instance.update_M(B);

      // M Step Update of eta
      eta = (a_eta - 1 + arma::accu(q_star) / 2) / (a_eta + b_eta - 2 + q * (q - 1) / 2);
      // M Step Update of Omega

      xi_star /= n; // QUIC works with re-scaled objective and we scale our penalty accordingly
      res_quic = cgquic(q, instance.S_Omega, instance.M, xi_star, eps, quic_max_iter);
      Omega = res_quic.slice(0);
      Sigma = res_quic.slice(1);
      instance.postprocessing(B,Sigma,Omega);
      
      if(verbose == 1) Rcout << "    updating residual matrices "<< endl;
      instance.update(mu_old,B,Sigma,Omega,n_rep,nskp);
      mu_old = instance.mu;

      // check convergence and whether we need to terminate early
      converged = 1;
      early_terminate = 0;
      omega_non_zero = 0;

      obj = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
      
      if (instance.s_eval(q - 1) / instance.s_eval(0) > s_max_condition)
        early_terminate = 1; // condition number of S is too large. Terminate early
      if ((obj - old_obj) / old_obj < eps)
        obj_counter++; // objective increased by less than 100*eps%, increment counter
      else
        obj_counter = 0; // if objective increases by more than 100*eps% re-set counter

      for (int j = 0; j < p; j++)
      {
        for (int k = 0; k < q; k++)
        {
          if ((B_old(j, k) == 0) & (B(j, k) != 0))
            converged = 0; // variable was added
          else if ((B_old(j, k) != 0) & (abs((B_old(j, k) - B(j, k)) / B_old(j, k)) > eps))
            converged = 0; // value of beta_{j,k} changed by more than 100*eps%
        }
      }
      for (int k = 0; k < q; k++)
      {
        for (int kk = k + 1; kk < q; kk++)
        {
          if ((Omega_old(k, kk) == 0) & (Omega(k, kk) != 0))
            converged = 0;
          else if ((Omega_old(k, kk) != 0) & (abs((Omega_old(k, kk) - Omega(k, kk)) / Omega_old(k, kk)) > eps))
            converged = 0;
          if (Omega(k, kk) != 0)
            omega_non_zero++;
        }
      }
      if (obj_counter == obj_counter_max)
        break;
      if (early_terminate == 1)
        break;
      if (converged == 1)
        break;
    } // closes while loop of main EM algorithm
    if (obj_counter == obj_counter_max)
      obj_term(current_index) = 1;
    early_term(current_index) = early_terminate;
    conv(current_index) = converged;
    if (verbose == 1)
    {
      if (early_terminate == 1)
        Rcout << "    Terminated early due to ill-conditioned S." << endl;
      if (obj_counter == obj_counter_max)
        Rcout << "    Terminated early due to marginal objective increase." << endl;
      if ((iter == max_iter) & (converged == 0))
        Rcout << "    Reached max iterations but failed to converge." << endl;
      Rcout << "    num B != 0 : " << arma::accu(B != 0) << " num Omega != 0 : " << omega_non_zero;
      Rcout << "Finished s = " << s << "t = " << t << endl;
    }
    // save to paths
    B_path.slice(current_index) = B;
    Omega_path.slice(current_index) = Omega;
    Sigma_path.slice(current_index) = Sigma;
    R_path.slice(current_index) = instance.R;
    tXR_path.slice(current_index) = instance.tXR;
    S_B_path.slice(current_index) = instance.S;
    S_Omega_path.slice(current_index) = instance.S_Omega;
    M_path.slice(current_index) = instance.M;
    theta_path(current_index) = theta;
    eta_path(current_index) = eta;
    obj_path(current_index) = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda_spike(L - 1), xi1, xi_spike(L - 1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
  }

  // Now we need to solve for s = 0
  s = 0;
  t = 0;
  for (t = 1; t < L; t++)
  {
    lambda0 = lambda_spike(s);
    xi0 = xi_spike(t);
    if (verbose == 1)
      Rcout << "Started s = " << s << " t = " << t << endl;
    current_index = s + L * t;
    down_index = s + L * (t - 1);
    if (early_term(down_index) == 1)
    {
      B = B_restart;
      Omega = Omega_restart;
      Sigma = Sigma_restart;
      instance = instance_restart;
      theta = theta_restart;
      eta = eta_restart;
    }
    else
    {
      B = B_path.slice(down_index);
      Omega = Omega_path.slice(down_index);
      Sigma = Sigma_path.slice(down_index);
      instance.R = R_path.slice(down_index);
      instance.tXR = tXR_path.slice(down_index);
      instance.S = S_B_path.slice(down_index);
      instance.S_Omega = S_Omega_path.slice(down_index);
      instance.M = M_path.slice(down_index);
      theta = theta_path(down_index);
      eta = eta_path(down_index);
    }
    prop_dir(current_index) = 2;
    iter = 0;
    converged = 0;
    early_terminate = 0;
    while (iter < max_iter)
    {
      iter++;
      if (iter % 50)
        Rcpp::checkUserInterrupt();
      B_old = B;
      Omega_old = Omega;
      old_obj = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);

      // E Step
      for (int k = 0; k < q; k++)
      {
        for (int kk = k + 1; kk < q; kk++)
        {
          tmp = 1.0 / (1.0 + (1.0 - eta) / eta * xi0 / xi1 * exp(-1.0 * abs(Omega(k, kk)) * (xi0 - xi1)));
          q_star(k, kk) = tmp;
          q_star(kk, k) = tmp;
        }
      }
      q_star.diag().zeros(); // to update eta, we need sum of q_star's easier to set diagonals equal to 0 and sum over entire matrix
      xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
      if (diag_penalty == 1)
        xi_star.diag().fill(xi1);
      else if (diag_penalty == 0)
        xi_star.diag().fill(0);
      
      // M Step Update of B and Theta
      
      update_B_theta(instance.n_B, p, q, B, instance.R, instance.tXR, instance.S, theta, Sigma, eta, instance.X, instance.tXX, lambda1, lambda0, xi1, xi0, diag_penalty, theta_hyper_params, eta_hyper_params, max_iter, eps, verbose);
      // update things related with B
      instance.update_M(B);
      // M Step Update of eta
      eta = (a_eta - 1 + arma::accu(q_star) / 2) / (a_eta + b_eta - 2 + q * (q - 1) / 2);
      // M Step Update of Omega
      xi_star /= n; // QUIC works with re-scaled objective and we scale our penalty accordingly

      res_quic = cgquic(q, instance.S_Omega, instance.M, xi_star, eps, quic_max_iter);
      Omega = res_quic.slice(0);
      Sigma = res_quic.slice(1);
      instance.postprocessing(B,Sigma,Omega);
      
      if(verbose == 1) Rcout << "    updating residual matrices "<< endl;
      instance.update(mu_old,B,Sigma,Omega,n_rep,nskp);
      mu_old = instance.mu;

      // check convergence and whether we need to terminate early
      converged = 1;
      early_terminate = 0;
      omega_non_zero = 0;

      obj = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
      
      if (instance.s_eval(q - 1) / instance.s_eval(0) > s_max_condition)
        early_terminate = 1; // condition number of S is too large. Terminate early
      if ((obj - old_obj) / old_obj < eps)
        obj_counter++; // objective increased by less than 100*eps%, increment counter
      else
        obj_counter = 0; // if objective increases by more than 100*eps% re-set counter

      for (int j = 0; j < p; j++)
      {
        for (int k = 0; k < q; k++)
        {
          if ((B_old(j, k) == 0) & (B(j, k) != 0))
            converged = 0; // variable was added
          else if ((B_old(j, k) != 0) & (abs((B_old(j, k) - B(j, k)) / B_old(j, k)) > eps))
            converged = 0; // value of beta_{j,k} changed by more than 100*eps%
        }
      }
      for (int k = 0; k < q; k++)
      {
        for (int kk = k + 1; kk < q; kk++)
        {
          if ((Omega_old(k, kk) == 0) & (Omega(k, kk) != 0))
            converged = 0;
          else if ((Omega_old(k, kk) != 0) & (abs((Omega_old(k, kk) - Omega(k, kk)) / Omega_old(k, kk)) > eps))
            converged = 0;
          if (Omega(k, kk) != 0)
            omega_non_zero++;
        }
      }
      if (obj_counter == obj_counter_max)
        break;
      if (early_terminate == 1)
        break;
      if (converged == 1)
        break;
    } // closes while loop of main EM algorithm
    if (obj_counter == obj_counter_max)
      obj_term(current_index) = 1;
    early_term(current_index) = early_terminate;
    conv(current_index) = converged;
    if (verbose == 1)
    {
      if (early_terminate == 1)
        Rcout << "    Terminated early due to ill-conditioned S." << endl;
      if (obj_counter == obj_counter_max)
        Rcout << "    Terminated early due to marginal objective increase." << endl;
      if ((iter == max_iter) & (converged == 0))
        Rcout << "    Reached max iterations but failed to converge." << endl;
      Rcout << "    num B != 0 : " << arma::accu(B != 0) << " num Omega != 0 : " << omega_non_zero;
      Rcout << "Finished s = " << s << "t = " << t << endl;
    }
    // save to paths
    B_path.slice(current_index) = B;
    Omega_path.slice(current_index) = Omega;
    Sigma_path.slice(current_index) = Sigma;
    R_path.slice(current_index) = instance.R;
    tXR_path.slice(current_index) = instance.tXR;
    S_B_path.slice(current_index) = instance.S;
    S_Omega_path.slice(current_index) = instance.S_Omega;
    M_path.slice(current_index) = instance.M;
    theta_path(current_index) = theta;
    eta_path(current_index) = eta;
    obj_path(current_index) = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda_spike(L - 1), xi1, xi_spike(L - 1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
  }
  // Now we're ready to solve with s > 0 and t > 0

  for (s = 1; s < L; s++)
  {
    for (t = 1; t < L; t++)
    {
      if (verbose == 1)
        Rcout << "Started s = " << s << " t = " << t << endl;
      lambda0 = lambda_spike(s);
      xi0 = xi_spike(t);
      current_index = s + L * t;

      left_index = s - 1 + L * t;
      down_index = s + L * (t - 1);
      diag_index = (s - 1) + L * (t - 1);

      if (early_term(left_index) == 1)
      { // solution to left was unstable. use the re-start
        B_left = B_restart;
        Omega_left = Omega_restart;
        Sigma_left = Sigma_restart;
        instance_left = instance_restart;
        theta_left = theta_restart;
        eta_left = eta_restart;
      }
      else
      {
        B_left = B_path.slice(left_index);
        Omega_left = Omega_path.slice(left_index);
        Sigma_left = Sigma_path.slice(left_index);
        instance_left.R = R_path.slice(left_index);
        instance_left.tXR = tXR_path.slice(left_index);
        instance_left.S = S_B_path.slice(left_index);
        instance_left.S_Omega = S_Omega_path.slice(left_index);
        instance_left.M = M_path.slice(left_index);
        theta_left = theta_path(left_index);
        eta_left = eta_path(left_index);
      }
      obj_left = cgobjective(n, p, q, instance_left.S_Omega, B_left, instance.X, Omega_left, Sigma_left, lambda1, lambda0, xi1, xi0, theta_left, eta_left, diag_penalty, theta_hyper_params, eta_hyper_params);
      if (early_term(down_index) == 1)
      {
        B_down = B_restart;
        Omega_down = Omega_restart;
        Sigma_down = Sigma_restart;
        instance_down = instance_restart;
        theta_down = theta_restart;
        eta_down = eta_restart;
      }
      else
      {
        B_down = B_path.slice(down_index);
        Omega_down = Omega_path.slice(down_index);
        Sigma_down = Sigma_path.slice(down_index);
        instance_down.R = R_path.slice(down_index);
        instance_down.tXR = tXR_path.slice(down_index);
        instance_down.S = S_B_path.slice(down_index);
        instance_down.S_Omega = S_Omega_path.slice(down_index);
        instance_down.M = M_path.slice(down_index);
        theta_down = theta_path(down_index);
        eta_down = eta_path(down_index);
      }
      obj_down = cgobjective(n, p, q, instance_down.S_Omega, B_down, instance.X, Omega_down, Sigma_down, lambda1, lambda0, xi1, xi0, theta_down, eta_down, diag_penalty, theta_hyper_params, eta_hyper_params);
      if (early_term(diag_index) == 1)
      {
        B_diag = B_restart;
        Omega_diag = Omega_restart;
        Sigma_diag = Sigma_restart;
        instance_diag = instance_restart;
        theta_diag = theta_restart;
        eta_diag = eta_restart;
      }
      else
      {
        B_diag = B_path.slice(diag_index);
        Omega_diag = Omega_path.slice(diag_index);
        Sigma_diag = Sigma_path.slice(diag_index);
        instance_diag.R = R_path.slice(diag_index);
        instance_diag.tXR = tXR_path.slice(diag_index);
        instance_diag.S = S_B_path.slice(diag_index);
        instance_diag.S_Omega = S_Omega_path.slice(diag_index);
        instance_diag.M = M_path.slice(diag_index);
        theta_diag = theta_path(diag_index);
        eta_diag = eta_path(diag_index);
      }
      obj_diag = cgobjective(n, p, q, instance_diag.S_Omega, B_diag, instance.X, Omega_diag, Sigma_diag, lambda1, lambda0, xi1, xi0, theta_diag, eta_diag, diag_penalty, theta_hyper_params, eta_hyper_params);

      if ((obj_left >= obj_down) & (obj_left >= obj_diag))
      {
        B = B_left;
        Omega = Omega_left;
        Sigma = Sigma_left;
        instance = instance_left;
        theta = theta_left;
        eta = eta_left;
        prop_dir(current_index) = 1;
        old_obj = obj_left;
      }
      else if ((obj_down > obj_left) & (obj_down >= obj_diag))
      {
        B = B_down;
        Omega = Omega_down;
        Sigma = Sigma_down;
        instance = instance_down;
        theta = theta_down;
        eta = eta_down;
        prop_dir(current_index) = 2;
        old_obj = obj_down;
      }
      else if ((obj_diag > obj_left) & (obj_diag > obj_down))
      {
        B = B_diag;
        Omega = Omega_diag;
        Sigma = Sigma_diag;
        instance = instance_diag;
        theta = theta_diag;
        eta = eta_diag;
        prop_dir(current_index) = 3;
        old_obj = obj_diag;
      }
      iter = 0;
      converged = 0;
      early_terminate = 0;
      while (iter < max_iter)
      {
        iter++;
        if (iter % 50)
          Rcpp::checkUserInterrupt();
        B_old = B;
        Omega_old = Omega;
        old_obj = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);

        // E Step
        for (int k = 0; k < q; k++)
        {
          for (int kk = k + 1; kk < q; kk++)
          {
            tmp = 1.0 / (1.0 + (1.0 - eta) / eta * xi0 / xi1 * exp(-1.0 * abs(Omega(k, kk)) * (xi0 - xi1)));
            q_star(k, kk) = tmp;
            q_star(kk, k) = tmp;
          }
        }
        q_star.diag().zeros(); // to update eta, we need sum of q_star's easier to set diagonals equal to 0 and sum over entire matrix
        xi_star = xi1 * q_star + xi0 * (1.0 - q_star);
        if (diag_penalty == 1)
          xi_star.diag().fill(xi1);
        else if (diag_penalty == 0)
          xi_star.diag().fill(0);

        
        // M Step Update of B and Theta
        update_B_theta(instance.n_B, p, q, B, instance.R, instance.tXR, instance.S, theta, Sigma, eta, instance.X, instance.tXX, lambda1, lambda0, xi1, xi0, diag_penalty, theta_hyper_params, eta_hyper_params, max_iter, eps, verbose);
        // update stuff related to B
        instance.update_M(B);
        // M Step Update of eta
        eta = (a_eta - 1 + arma::accu(q_star) / 2) / (a_eta + b_eta - 2 + q * (q - 1) / 2);
        // M Step Update of Omega
        xi_star /= n; // QUIC works with re-scaled objective and we scale our penalty accordingly
        res_quic = cgquic(q, instance.S_Omega, instance.M, xi_star, eps, quic_max_iter);
        Omega = res_quic.slice(0);
        Sigma = res_quic.slice(1);
        instance.postprocessing(B,Sigma,Omega);
        
        if(verbose == 1) Rcout << "    updating residual matrices "<< endl;
        instance.update(mu_old,B,Sigma,Omega,n_rep,nskp);
        mu_old = instance.mu;

        // check convergence and whether we need to terminate early
        converged = 1;
        early_terminate = 0;
        omega_non_zero = 0;

        obj = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda0, xi1, xi0, theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
        
        if (instance.s_eval(q - 1) / instance.s_eval(0) > s_max_condition)
          early_terminate = 1; // condition number of S is too large. Terminate early
        if ((obj - old_obj) / old_obj < eps)
          obj_counter++; // objective increased by less than 100*eps%, increment counter
        else
          obj_counter = 0; // if objective increases by more than 100*eps% re-set counter

        for (int j = 0; j < p; j++)
        {
          for (int k = 0; k < q; k++)
          {
            if ((B_old(j, k) == 0) & (B(j, k) != 0))
              converged = 0; // variable was added
            else if ((B_old(j, k) != 0) & (abs((B_old(j, k) - B(j, k)) / B_old(j, k)) > eps))
              converged = 0; // value of beta_{j,k} changed by more than 100*eps%
          }
        }
        for (int k = 0; k < q; k++)
        {
          for (int kk = k + 1; kk < q; kk++)
          {
            if ((Omega_old(k, kk) == 0) & (Omega(k, kk) != 0))
              converged = 0;
            else if ((Omega_old(k, kk) != 0) & (abs((Omega_old(k, kk) - Omega(k, kk)) / Omega_old(k, kk)) > eps))
              converged = 0;
            if (Omega(k, kk) != 0)
              omega_non_zero++;
          }
        }
        if (obj_counter == obj_counter_max)
          break;
        if (early_terminate == 1)
          break;
        if (converged == 1)
          break;
      } // closes while loop of main EM algorithm
      if (obj_counter == obj_counter_max)
        obj_term(current_index) = 1;
      early_term(current_index) = early_terminate;
      conv(current_index) = converged;
      if (verbose == 1)
      {
        if (early_terminate == 1)
          Rcout << "    Terminated early due to ill-conditioned S." << endl;
        if (obj_counter == obj_counter_max)
          Rcout << "    Terminated early due to marginal objective increase." << endl;
        if ((iter == max_iter) & (converged == 0))
          Rcout << "    Reached max iterations but failed to converge." << endl;
        Rcout << "    num B != 0 : " << arma::accu(B != 0) << " num Omega != 0 : " << omega_non_zero;
        Rcout << "Finished s = " << s << "t = " << t << endl;
      }
      // save to paths
      B_path.slice(current_index) = B;
      Omega_path.slice(current_index) = Omega;
      Sigma_path.slice(current_index) = Sigma;
      R_path.slice(current_index) = instance.R;
      tXR_path.slice(current_index) = instance.tXR;
      S_B_path.slice(current_index) = instance.S;
      S_Omega_path.slice(current_index) = instance.S_Omega;
      M_path.slice(current_index) = instance.M;
      theta_path(current_index) = theta;
      eta_path(current_index) = eta;
      //obj_path(current_index) = cgobjective(n, p, q, S_Omega, B, Omega, Sigma,lambda1, lambda_spike(L-1), xi1, xi_spike(L-1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params);
      obj_path(current_index) = cgobjective(n, p, q, instance.S_Omega, B, instance.X, Omega, Sigma, lambda1, lambda_spike(L - 1), xi1, xi_spike(L - 1), theta, eta, diag_penalty, theta_hyper_params, eta_hyper_params); // compute at the maximal value
    }                                                                                                                                                                                                // closes the loop over t
  }                                                                                                                                                                                                  // closes the loop over s

  int time_end = time(&tp);
  List results;
  for (int s = 0; s < L; s++)
  {
    for (int t = 0; t < L; t++)
    {
      current_index = s + L * t;
      tmp_B = B_path.slice(current_index);
      tmp_B.each_col() /= x_col_weights;
      B_path.slice(current_index) = tmp_B;
      alpha_path.col(current_index) = instance.mu - tmp_B.t() * mu_x; // compute the intercept
    }
  }
  alpha = instance.mu - tmp_B.t() * mu_x;
  results["alpha"] = alpha;
  results["B"] = B_path.slice(L * L - 1);
  results["Omega"] = Omega_path.slice(L * L - 1);
  results["Sigma"] = Sigma_path.slice(L * L - 1);
  //results["R"] = R_path.slice(L*L-1);
  //results["S"] = S_path.slice(L*L-1);
  results["theta"] = theta_path(L * L - 1);
  results["eta"] = eta_path(L * L - 1);
  results["alpha_path"] = alpha_path;
  results["B_path"] = B_path;
  results["Omega_path"] = Omega_path;
  results["Sigma"] = Sigma_path;
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
  if(verbose == 1) Rcout << "done" << endl;
  return (results);
}

}


#endif