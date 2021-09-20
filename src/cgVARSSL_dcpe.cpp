
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
List cgVARSSL_dcpe(
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
  int p = q;
    
  arma::mat X = Y;
  X.shed_row(n-1);

  // Center columns of X and Y
  // re-scale the re-centered columns of X to have norm sqrt(n)

  double tmp_mu_x = 0.0;
  double tmp_weight_x = 0.0;
  double tmp_mu_y = 0.0;
  arma::vec x_col_weights(p);
  arma::vec mu_x(p);
  arma::vec mu_y(q);
  for (int j = 0; j < p; j++)
  {
    tmp_mu_x = mean(X.col(j));
    X.col(j) -= tmp_mu_x;
    tmp_weight_x = norm(X.col(j)) / sqrt(n);
    X.col(j) /= tmp_weight_x;
    mu_x(j) = tmp_mu_x;
    x_col_weights(j) = tmp_weight_x;
  }
  for (int k = 0; k < q; k++)
  {
    tmp_mu_y = mean(Y.col(k));
    Y.col(k) -= tmp_mu_y;
    mu_y(k) = tmp_mu_y;
  }
  int quic_max_iter = 5 * max_iter;
  arma::mat Y_init = Y.row(0); // the first time step, not needed in estimating B
  Y.shed_row(0);// remove it 
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

  // initialize some useful matrices
  arma::mat tYY = Y.t() * Y + Y_init.t() * Y_init;
  arma::mat tXX = X.t() * X;
  arma::mat S_Omega = tYY / n; // S used for cgQUIC

  // initialize our parameters
  arma::vec alpha(q); // will hold the intercept. only computed at the end of all of the loops
  arma::mat B(p, q);
  B.zeros();
  arma::mat B_old = B; // used internally
  arma::mat tmp_B = B; // we will re-scale the working B matrix when we save the path
  arma::mat Omega(q, q);
  Omega.eye();
  //Omega = inv_sympd(S_Omega);
  arma::mat Omega_old = Omega; // used internally to check convergence in the EM algorithm
  arma::mat Sigma(q, q);
  Sigma.eye();
  //Sigma = S_Omega;
  arma::mat q_star(q, q);
  q_star.fill(0.5);
  arma::mat xi_star(q, q);
  xi_star.fill(xi1);
  // June 27: penalize diagonals
  xi_star.diag().fill(xi1);
  xi_star /= n; // the re-scaling is for compatability with QUIC.
  //xi_star.diag().zeros(); // for now, don't penalize the diagonals
  double tmp;             //holds value of q_star
  int omega_non_zero = 0; // counter for the number of non-zero off-diagonal elements of Omega's

  arma::mat mu_mat = X * B;
  arma::mat R = Y * Omega - mu_mat;
  arma::mat tXR = X.t() * R;
  arma::mat tRR = R.t() * R;
  arma::mat S_B = tRR / (n-1); // the copy of S used for updating B
  arma::mat M = mu_mat.t() * mu_mat / n;

  if (verbose == 1)
    Rcout << "Initialized R, tXR, S" << endl;

  // if we ever fail to converge, we should re-set the values of B and Omega
  arma::mat B_reset = B;
  arma::mat Omega_reset = Omega;
  arma::mat Sigma_reset = Sigma;
  arma::mat R_reset = R;
  arma::mat tXR_reset = tXR;
  arma::mat tRR_reset = tRR;
  arma::mat S_B_reset = S_B;
  arma::mat S_Omega_reset = S_Omega;
  arma::mat M_reset = M;
  double theta_reset = theta;
  double eta_reset = eta;

  //Rcout << "preparing QUIC working parameters" << endl;
  // initialize stuff for QUIC
  cube res_quic(q, q, 2); // cube to hold values from QUIC. 1st one is for Omega, 2nd is for Sigma
  //double* S_ptr = S.memptr();
  //double* Xi_ptr = xi_star.memptr();

  double a_eta = eta_hyper_params(0);
  double b_eta = eta_hyper_params(1);

  int converged = 1;
  int iter = 0;

  // Store the conditional modes
  cube B0_path(p, q, L);
  cube Omega0_path(q, q, L);
  arma::vec theta0_path(L);
  arma::vec eta0_path(L);
  if (verbose == 1)
    Rcout << "Starting Step 1: Estimating conditional modes of B" << endl;
  time_t tp;
  int time_start = time(&tp);
  for (int l = 0; l < L; l++)
  {
    lambda0 = lambda_spike(l);
    if (verbose == 1)
      Rcout << setprecision(6) << "Starting lambda0 = " << lambda0 << " num B non-zero = " << accu(B != 0) << " theta = " << theta << endl;
    update_B_theta(n-1, p, q, B, R, tXR, S_B, theta, Sigma, eta, X, tXX, lambda1, lambda0, xi1, xi0, diag_penalty, theta_hyper_params, eta_hyper_params, max_iter, eps, verbose);

    tmp_B = B;
    tmp_B.each_col() / x_col_weights;
    B0_path.slice(l) = tmp_B;
    theta0_path(l) = theta;
  }
  mu_mat = X * B;
  M = mu_mat.t() * mu_mat / n;
  // for debugging
  //Rcout << "checking S_B" <<endl;
  //Rcout << S_B << endl;
  //R = Y * Omega - mu_mat;
  //tXR = X.t() * R;
  //tRR = R.t() * R;
  //S_B = tRR/n;
  //Rcout << S_B << endl;
  //Rcout << "M" << endl;
  //Rcout << M << endl;
  // end testing

  if (verbose == 1)
    Rcout << "starting step 2 now" << endl;
  for (int l = 0; l < L; l++)
  {
    xi0 = xi_spike(l);
    omega_non_zero = 0;
    for (int k = 0; k < q; k++)
    {
      for (int kk = k + 1; kk < q; kk++)
      {
        if (Omega(k, kk) != 0)
          omega_non_zero++;
      }
    }
    if (verbose == 1)
      Rcout << "Starting xi0 = " << xi0 << " num Omega non-zero = " << omega_non_zero << " eta = " << eta << endl;
    iter = 0;
    while (iter < max_iter)
    {
      iter++;
      Omega_old = Omega;
      // E-step
      if (xi0 == xi1)
      {
        q_star.fill(1.0);
      }
      else
      {
        for (int k = 0; k < q; k++)
        {
          for (int kk = k + 1; kk < q; kk++)
          {
            tmp = 1.0 / (1.0 + (1.0 - eta) / eta * xi0 / xi1 * exp(-1.0 * abs(Omega(k, kk)) * (xi0 - xi1)));
            q_star(k, kk) = tmp;
            q_star(kk, k) = tmp;
          }
        }
      }
      q_star.diag().zeros();
      xi_star = xi1 * q_star + xi0 * (1 - q_star);
      if (diag_penalty == 1)
        xi_star.diag().fill(xi1);
      else if (diag_penalty == 0)
        xi_star.diag().fill(0);
      // M-step update of eta
      eta = (a_eta - 1 + accu(q_star) / 2) / (a_eta + b_eta - 2 + q * (q - 1) / 2);
      // M-step update of Omega
      xi_star /= n; // QUIC needs everything to be scaled
      //Rcout << "xi_star" << endl;
      //Rcout << xi_star << endl;
      res_quic = cgquic(q, S_Omega, M, xi_star, eps, quic_max_iter);
      Omega = res_quic.slice(0);
      Sigma = res_quic.slice(1);

      // check convergence
      converged = 1;
      omega_non_zero = 0;
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
      if (converged == 1)
      {
        //Rcout << "step 2 Omega converged" << endl;
        break;
      }
    } // closes while loop (i.e. the main EM)
    if ((iter == max_iter) & (converged == 0))
    {
      Rcout << "Omega did not converge. Re-setting the values" << endl;
      Omega = Omega_reset;
      Sigma = Sigma_reset;
      eta = eta_reset;
    }
    else if (converged == 1)
    {
      Omega_reset = Omega;
      Sigma_reset = Sigma;
      eta_reset = eta;
      //Rcout << "flag" <<endl;
    }
    Omega0_path.slice(l) = Omega;

  } // closes for loop over l
  // updating things related with Omega
  R = Y * Omega - mu_mat;
  tXR = X.t() * R;
  tRR = R.t() * R;
  S_B = tRR / (n-1); // the copy of S used for updating B
  //Rcout << Omega << endl;
  if (verbose == 1)
    Rcout << "starting Step 3" << endl;
  B_reset = B;
  Omega_reset = Omega;
  Sigma_reset = Sigma;
  R_reset = R;
  tXR_reset = tXR;
  S_B_reset = S_B;
  S_Omega_reset = S_Omega;
  M_reset = M;
  theta_reset = theta;
  eta_reset = eta;

  lambda0 = lambda_spike(L - 1);
  xi0 = xi_spike(L - 1);
  iter = 0;
  converged = 1;
  while (iter < max_iter)
  {
    iter++;
    B_old = B;
    Omega_old = Omega;

    update_B_theta(n-1, p, q, B, R, tXR, S_B, theta, Sigma, eta, X, tXX, lambda1, lambda0, xi1, xi0, diag_penalty, theta_hyper_params, eta_hyper_params, max_iter, eps, verbose);

    mu_mat = X * B;
    M = mu_mat.t() * mu_mat / n;
    // for debugging
    //Rcout << "checking S_B" <<endl;
    //Rcout << S_B << endl;
    //R = Y * Omega - mu_mat;
    //tXR = X.t() * R;
    //tRR = R.t() * R;
    //S_B = tRR/n;
    //Rcout << S_B << endl;
    //Rcout << "M" << endl;
    //Rcout << M <<endl;
    // end testing

    // now update Omega and eta
    for (int k = 0; k < q; k++)
    {
      for (int kk = k + 1; kk < q; kk++)
      {
        tmp = 1.0 / (1.0 + (1.0 - eta) / eta * xi0 / xi1 * exp(-1.0 * abs(Omega(k, kk)) * (xi0 - xi1)));
        q_star(k, kk) = tmp;
        q_star(kk, k) = tmp;
      }
    }
    q_star.diag().zeros();
    xi_star = xi1 * q_star + xi0 * (1 - q_star);
    if (diag_penalty == 1)
      xi_star.diag().fill(xi1);
    else if (diag_penalty == 0)
      xi_star.diag().fill(0);
    // M-step update of eta
    eta = (a_eta - 1 + accu(q_star) / 2) / (a_eta + b_eta - 2 + q * (q - 1) / 2);
    // M-step update of Omega
    xi_star /= n;
    res_quic = cgquic(q, S_Omega, M, xi_star, eps, quic_max_iter);
    Omega = res_quic.slice(0);
    Sigma = res_quic.slice(1);

    // updating things related with Omega
    R = Y * Omega - mu_mat;
    tXR = X.t() * R;
    tRR = R.t() * R;
    S_B = tRR / (n-1); // the copy of S used for updating B

    converged = 1;
    for (int j = 0; j < p; j++)
    {
      for (int k = 0; k < q; k++)
      {
        if ((B_old(j, k) == 0) & (B(j, k) != 0))
          converged = 0;
        else if ((B_old(j, k) != 0) & (abs((B_old(j, k) - B(j, k)) / B_old(j, k)) > eps))
          converged = 2;
      }
    }
    omega_non_zero = 0;
    for (int k = 0; k < q; k++)
    {
      for (int kk = k + 1; kk < q; kk++)
      {
        if ((Omega_old(k, kk) == 0) & (Omega(k, kk) != 0))
          converged = 3;
        else if ((Omega_old(k, kk) != 0) & (abs((Omega_old(k, kk) - Omega(k, kk)) / Omega_old(k, kk)) > eps))
        {
          converged = 4;
          //Rcout << "k" << k << " kk " << kk << endl;
        }
        if (Omega(k, kk) != 0)
          omega_non_zero++;
      }
    }
    if (verbose == 1)
      Rcout << "iter = " << iter << " num nonzero B: " << accu(B != 0) << " num nonzero Omega: " << omega_non_zero << " theta = " << theta << " eta = " << eta << endl;
    if (converged == 1)
      break;
  }
  if ((iter == max_iter) & (converged != 1))
  {
    // did not converge!
    Rcout << "Omega and B did not converge with code " << converged << endl;
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
