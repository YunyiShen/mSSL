#ifndef WORKINGPARAM_MPARAM_H
#define WORKINGPARAM_MPARAM_H

// [[Rcpp::depends(RcppArmadillo)]]

#include "wparam.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

namespace workingparam{

class mWorkingParam : public WorkingParam{
    public:
        mWorkingParam() = default;
        mWorkingParam(arma::mat X, arma::mat Y) : WorkingParam(X,Y){}
        inline void update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        int n_rep, int nskp = 5);
};

inline void mWorkingParam::update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,                    
                        int n_rep, int nskp){

    arma::mat XB = X * B_t;
    XB.each_row() += mu_t.t();
    R = Y-XB;
    tRR = R.t() * R;
    S = tRR/n_B;
    tXR = X.t() * R;

}



}
#endif