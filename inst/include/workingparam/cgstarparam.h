#ifndef WORKINGPARAM_CGSTARPARAM_H
#define WORKINGPARAM_CGSTARPARAM_H


// [[Rcpp::depends(RcppArmadillo)]]
#include <linconGaussR.h> 
#include "wparam.h"

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace linconGaussR;

namespace workingparam{

class cgstarWorkingParam: public cgWorkingParam {
    
    public: 
        arma::mat lower;
        arma::mat upper;
        cgstarWorkingParam() = default;
        cgstarWorkingParam(arma::mat X, arma::mat lower, arma::mat upper): cgWorkingParam(X,lower),lower(lower),upper(upper){}

        inline void update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        const arma::mat &Omega_t,
                        int n_rep, int nskp = 5);
};

inline void cgstarWorkingParam::update(
                        const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        const arma::mat &Omega_t,
                        int n_rep, int nskp){
    int n = lower.n_rows;
    int q = upper.n_cols;
    arma::mat XB = X * B_t;
    XB.each_row() += mu_t.t();
    M = XB.t() * XB / n;
    arma::vec x_init(q);
    arma::mat C = arma::chol(Sigma_t);
    S.zeros(q,q);
    S_Omega.zeros(q,q);
    mu.zeros(q);
    R.zeros(n,q);
    tXR.zeros(q,q);
    s_eval.zeros(q);

    // several helpers
    arma::rowvec meani;
    arma::uvec finites;
    arma::vec x_init_temp;

    for(int i = 0 ; i < n ; i++){
        meani = XB.row(i) * Sigma_t;
        Rcpp::checkUserInterrupt();
        arma::vec b(2*q);
        arma::mat A(2*q,q);
        b.rows(0,q-1) = trans(meani-lower.row(i));
        b.rows(q,2*q-1) = trans(upper.row(i)-meani);
        finites = find_finite(b);
        A.rows(0,q-1) = C.t();
        A.rows(q,2*q-1) = -C.t();
        b = b.rows(finites);
        A = A.rows(finites);
        x_init_temp = trans( lower.row(i) );
        x_init = trans(upper.row(i));
        x_init_temp(find_nonfinite(x_init_temp)) = x_init(find_nonfinite(x_init_temp))-1;
        x_init(find_nonfinite(x_init)) = x_init_temp(find_nonfinite(x_init))+1;
        x_init += x_init_temp;
        x_init *= 0.5;
        x_init = arma::solve(C.t(), x_init-arma::trans(meani));// this serves as the initial point
        LinearConstraints lincon(A,b,true);
        //Rcout << "sampling " << i << "th sample" << endl;
        EllipticalSliceSampler sampler(n_rep + 1,lincon,nskp,x_init);
        sampler.run();
        //Rcout << "      end" << endl;
        arma::mat resi = sampler.loop_state.samples;
        resi.shed_row(0);// we will remove the initial points
        resi = resi * C;
        // the part you don't want to transform Y*:
        resi.each_row() += meani;
        S_Omega += resi.t() * resi;

        // now we have to transform it:
        resi = resi * Omega_t;
        resi.each_row() -= XB.row(i);
        R.row(i) = mean(resi);
        S += resi.t() * resi;
        resi.each_row() += mu_t.t();
        mu += arma::trans(arma::mean(resi)); 
    }

    S /= (n*n_rep);
    S_Omega /= (n*n_rep);
    mu /= n;
    R.each_row() += mu_t.t();
    R.each_row() -= mu.t();
    tXR = X.t() * R;
    tRR = S*n;
    
    s_eval = eig_sym(S);

}


}

#endif