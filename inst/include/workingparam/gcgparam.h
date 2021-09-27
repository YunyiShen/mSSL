#ifndef WORKINGPARAM_CGPARAM_H
#define WORKINGPARAM_CGPARAM_H

// [[Rcpp::depends(RcppArmadillo)]]

#include "wparam.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

namespace workingparam{

class gcgWorkingParam : public cgWorkingParam{
    public:
        gcgWorkingParam() = default;
        gcgWorkingParam(arma::mat X, arma::mat Y) : cgWorkingParam(X,Y){}
        inline void update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        const arma::mat &Omega_t,
                        int n_rep, int nskp = 5);
};

inline void gcgWorkingParam::update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        const arma::mat &Omega_t,
                        int n_rep, int nskp){

    arma::mat YOmega = Y * Omega_t;
    arma::mat XB = X * B_t;
    XB.each_row() += mu_t.t();
    R = YOmega-XB;
    tRR = R.t() * R;
    S = tRR/n_B;
    tXR = X.t() * R;
    mu = trans(mean(YOmega));

}



}
#endif