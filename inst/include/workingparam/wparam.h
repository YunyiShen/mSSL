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
        int n_B;
        int n_Omega; 
        WorkingParam() = default;
        WorkingParam(arma::mat SS, arma::mat RR,arma::mat tXRR,arma::mat tXXx,arma::vec muu,int n){
            S = SS;
            tRR = SS*n;
            R = RR;
            tXR = tXRR;
            tXX = tXXx;
            mu = muu;
            s_eval = eig_sym(S);
            n_B = R.n_rows;
            n_Omega = n_B;
        }
};

class cgWorkingParam{
    public:
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
        cgWorkingParam(arma::mat SS, arma::mat RR,arma::mat tXRR,arma::mat tXXx,arma::vec muu,int n){
            S = SS;
            tRR = SS*n;
            R = RR;
            tXR = tXRR;
            tXX = tXXx;
            mu = muu;
            S_Omega = SS;
            M = SS;
            s_eval = eig_sym(S);
            n_B = R.n_rows;
            n_Omega = n_B;
        }
        inline void update_M(const arma::mat &X, const arma::mat &B_new);
};

inline void cgWorkingParam::update_M(const arma::mat &X,const arma::mat &B_new){
    int n = X.n_rows;
    arma::mat XB = X * B_new;
    XB.each_row() += mu.t();
    M = XB.t() * XB / n;
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