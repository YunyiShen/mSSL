#ifndef WORKINGPARAM_CGVARPARAM_H
#define WORKINGPARAM_CGVARPARAM_H

// [[Rcpp::depends(RcppArmadillo)]]

#include "wparam.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

namespace workingparam{

class cgVARWorkingParam : public cgWorkingParam{
    public:
        arma::mat Y_B; // the sheded Y
        cgVARWorkingParam() = default;
        cgVARWorkingParam(arma::mat Y1){
            Y_B = Y1;
            Y_B.shed_row(0);
            X = Y1;
            X.shed_row(X.n_rows-1);
            Y = Y1;
            n_B = X.n_rows;
            n_Omega = Y.n_rows;
            tRR = Y_B.t() * Y_B;
            S = tRR/n_B;
            R = Y_B;
            tXR = X.t() * R;
            tXX = X.t() * X;
            mu = trans(mean(Y));
            S_Omega = Y.t()*Y/n_Omega;
            M = S;
            s_eval = eig_sym(S);
        }
        inline void update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        const arma::mat &Omega_t,
                        int n_rep, int nskp = 5);
};

inline void cgVARWorkingParam::update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        const arma::mat &Omega_t,
                        int n_rep, int nskp){

    //Rcout << "flag" << endl;
    arma::mat YOmega = Y_B * Omega_t;
    arma::mat XB = X * B_t;
    XB.each_row() += mu_t.t();
    R = YOmega-XB;
    //Rcout << "flag2" << endl;
    tRR = R.t() * R;
    S = tRR/n_B;
    tXR = X.t() * R;
    mu = trans(mean(Y)*Omega_t);

}



}
#endif