#include <workingparam.h>

//' scale covariance and precision in place
//' This routine will scale the covariance matrix to have unit diagonal in place and also properly scale the precision.
//' @param Sigma a matrix that is the covariance matrix, will be scaled to have unit diagonal
//' @param Omega should be the inverse of Sigma
//' @return no return, scaling will be done in place
// [[Rcpp::export]]
void unitdiag2(arma::mat & Sigma, arma::mat &Omega){
    workingparam::unitdiag(Sigma, Omega);
    return;
}