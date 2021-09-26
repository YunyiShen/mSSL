
#include <mSSL.h>

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace mSSL;
using namespace quic;
using namespace workingparam;

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
  mWorkingParam Worker(X,Y);
  List results = mSSL::mSSL_dpe<mWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                s_max_condition,obj_counter_max,verbose);
  return results;
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
  mWorkingParam Worker(X,Y);
  List results = mSSL::mSSL_dcpe<mWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                verbose);
  return results;
}



// [[Rcpp::export]]
List fstarSSL_dpe(arma::mat X,
              arma::mat lower,
              arma::mat upper,
              List lambdas,
              List xis,
              arma::vec theta_hyper_params,
              arma::vec eta_hyper_params,
              int diag_penalty,
              int max_iter,
              double eps,
              int s_max_condition,
              int obj_counter_max,
              int verbose, int nrep=200, int nskp=1)
{
  starWorkingParam Worker(X,lower,upper);
  List results = mSSL::mSSL_dpe<starWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                s_max_condition,obj_counter_max,verbose, nrep, nskp);
  return results;
}

// [[Rcpp::export]]
List fstarSSL_dcpe(arma::mat X,
              arma::mat lower,
              arma::mat upper,
              List lambdas,
              List xis,
              arma::vec theta_hyper_params,
              arma::vec eta_hyper_params,
              int diag_penalty,
              int max_iter,
              double eps,
              int verbose, int nrep = 200, int nskp = 1)
{
  starWorkingParam Worker(X,lower, upper);
  List results = mSSL::mSSL_dcpe<starWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                verbose, nrep, nskp);
  return results;
}

// [[Rcpp::export]]
List mpSSL_dpe(arma::mat X,
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
              int verbose, int nrep=200, int nskp=1)
{
  probitWorkingParam Worker(X,Y);
  List results = mSSL::mSSL_dpe<probitWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                s_max_condition,obj_counter_max,verbose, nrep, nskp);
  return results;
}


// [[Rcpp::export]]
List mpSSL_dcpe(arma::mat X,
              arma::mat Y,
              List lambdas,
              List xis,
              arma::vec theta_hyper_params,
              arma::vec eta_hyper_params,
              int diag_penalty,
              int max_iter,
              double eps,
              int verbose, int nrep, int nskp)
{
  probitWorkingParam Worker(X,Y);
  List results = mSSL::mSSL_dcpe<probitWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                verbose, nrep, nskp);
  return results;
}



