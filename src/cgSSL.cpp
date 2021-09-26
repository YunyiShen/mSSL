
#include <mSSL.h>

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace mSSL;
using namespace quic;
using namespace workingparam;

// [[Rcpp::export]]
List cgSSL_dpe(arma::mat X,
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
  gcgWorkingParam Worker(X,Y);
  List results = mSSL::cgSSL_dpe<gcgWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                s_max_condition,obj_counter_max,verbose);
  return results;
}


// [[Rcpp::export]]
List cgSSL_dcpe(arma::mat X,
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
  gcgWorkingParam Worker(X,Y);
  List results = mSSL::cgSSL_dcpe<gcgWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                verbose);
  return results;
}



// [[Rcpp::export]]
List cgfstarSSL_dpe(arma::mat X,
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
  cgstarWorkingParam Worker(X,lower,upper);
  List results = mSSL::cgSSL_dpe<cgstarWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                s_max_condition,obj_counter_max,verbose, nrep, nskp);
  return results;
}

// [[Rcpp::export]]
List cgfstarSSL_dcpe(arma::mat X,
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
  cgstarWorkingParam Worker(X,lower, upper);
  List results = mSSL::cgSSL_dcpe<cgstarWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                verbose, nrep, nskp);
  return results;
}

// [[Rcpp::export]]
List cgpSSL_dpe(arma::mat X,
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
  cgprobitWorkingParam Worker(X,Y);
  List results = mSSL::cgSSL_dpe<cgprobitWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                s_max_condition,obj_counter_max,verbose, nrep, nskp);
  return results;
}


// [[Rcpp::export]]
List cgpSSL_dcpe(arma::mat X,
              arma::mat Y,
              List lambdas,
              List xis,
              arma::vec theta_hyper_params,
              arma::vec eta_hyper_params,
              int diag_penalty,
              int max_iter,
              double eps,
              int verbose, int nrep = 200, int nskp = 1)
{
  cgprobitWorkingParam Worker(X,Y);
  List results = mSSL::cgSSL_dcpe<cgprobitWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                verbose, nrep, nskp);
  return results;
}

// [[Rcpp::export]]
List cgVARSSL_dpe(arma::mat Y,
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
  cgVARWorkingParam Worker(Y);
  List results = mSSL::cgSSL_dpe<cgVARWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                s_max_condition,obj_counter_max,verbose);
  return results;
}


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
  cgVARWorkingParam Worker(Y);
  List results = mSSL::cgSSL_dcpe<cgVARWorkingParam>(Worker, lambdas, 
                xis,theta_hyper_params,eta_hyper_params,
                diag_penalty,max_iter,eps,
                verbose);
  return results;
}


