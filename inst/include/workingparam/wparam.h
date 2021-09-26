#ifndef WORKINGPARAM_WPARAM_H
#define WORKINGPARAM_WPARAM_H

#include <linconGaussR.h> 

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace linconGaussR;

namespace workingparam{

class WorkingParam{
    public:
        arma::mat X;
        arma::mat Y;
        arma::mat S;
        arma::mat R;
        arma::mat tXX;
        arma::mat tXR;
        arma::mat tRR;
        arma::vec mu;
        arma::vec s_eval;
        int n_B;
        int n_Omega; 
        WorkingParam() = default;
        WorkingParam(arma::mat X, arma::mat Y): X(X),Y(Y){
            n_B = X.n_rows;
            n_Omega = Y.n_rows;
            tRR = Y.t() * Y;
            S = tRR/n_B;
            R = Y;
            tXR = X.t() * R;
            tXX = X.t() * X;
            mu = trans(mean(Y));
            s_eval = eig_sym(S);
        }
        inline void update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        int n_rep, int nskp){
            Rcpp::stop("Update method not implemented! Implement the residual matrix updating method");
        }
        inline void postprocessing(
                        arma::mat &B,
                        arma::mat &Sigma,
                        arma::mat &Omega){
            return;
        }
};

class cgWorkingParam{
    public:
        arma::mat X;
        arma::mat Y;
        arma::mat S;
        arma::mat R;
        arma::mat tXX;
        arma::mat tXR;
        arma::mat tRR;
        arma::mat S_Omega;
        arma::mat M;
        arma::vec mu;
        arma::vec s_eval;
        int n_B;
        int n_Omega;
        cgWorkingParam() = default;
        cgWorkingParam(arma::mat X, arma::mat Y):X(X),Y(Y){
            n_B = X.n_rows;
            n_Omega = Y.n_rows;
            tRR = Y.t() * Y;
            S = tRR/n_B;
            R = Y;
            tXR = X.t() * R;
            tXX = X.t() * X;
            mu = trans(mean(Y));
            S_Omega = S;
            M = S;
            s_eval = eig_sym(S);
        }
        inline void update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        int n_rep, int nskp){
            Rcpp::stop("Update method not implemented! Implement the residual matrix updating method");
        }
        inline void update_M(const arma::mat &B_new);
        inline void postprocessing(
                        arma::mat &B,
                        arma::mat &Sigma,
                        arma::mat &Omega){
            return;
        }
};

inline void cgWorkingParam::update_M(const arma::mat &B_new){
    arma::mat XB = X * B_new;
    XB.each_row() += mu.t();
    M = XB.t() * XB / n_B;
}

inline void unitdiag(arma::mat & Sigma, arma::mat &Omega){
    arma::vec scaling = Sigma.diag();
    scaling = arma::sqrt(scaling);
    Omega.each_col() %= scaling;
    Omega.each_row() %= scaling.t();
    scaling = 1/scaling;
    Sigma.each_row() %= scaling.t();
    Sigma.each_col() %= scaling;
    return;
}


}




#endif