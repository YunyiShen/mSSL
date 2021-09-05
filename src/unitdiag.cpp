// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

//[[Rcpp::export]]
void unitdiag(arma::mat & Sigma, arma::mat &Omega){
    arma::vec scaling = Sigma.diag();
    scaling = arma::sqrt(scaling);
    Omega.each_col() %= scaling;
    Omega.each_row() %= scaling.t();
    scaling = 1/scaling;
    Sigma.each_row() %= scaling.t();
    Sigma.each_col() %= scaling;
    return;
}
