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
        arma::mat S;
        arma::mat R;
        arma::mat tXX;
        arma::mat tXR;
        arma::mat tRR;
        arma::vec mu;
        arma::vec s_eval;
        WorkingParam() = default;
        WorkingParam(arma::mat SS, arma::mat RR,arma::mat tXRR,arma::mat tXXx,arma::vec muu,int n){
            S = SS;
            tRR = SS*n;
            R = RR;
            tXR = tXRR;
            tXX = tXXx;
            mu = muu;
            s_eval = eig_sym(S);
        }
    };

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