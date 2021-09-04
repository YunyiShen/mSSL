// gSSL_obj.h
// [[Rcpp:depends(RcppArmadillo)]]

#ifndef GSSL_OBJ_H
#define GSSL_OBJ_H

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
#endif

double gSSL_objective(int n, int q, arma::mat S, arma::mat Omega, double xi1, double xi0, double eta, int diag_penalty, arma::vec eta_hyper_params){
  
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

