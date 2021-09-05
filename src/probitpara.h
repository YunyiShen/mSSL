#ifndef PROBITPARA_H
#define PROBITPARA_H


// [[Rcpp::depends(RcppArmadillo)]]
#include "active_intersections.h"
#include "angle_sampler.h"
#include "ellipse.h"
#include "linear_constraints.h"
#include "loop.h"
#include "elliptical_slice_sampling.h"

#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;



class probitWorkingParam{
    public:
        arma::mat S;
        arma::mat R;
        arma::mat tXX;
        arma::mat tXR;
        arma::mat tRR;
        arma::vec mu;
        arma::vec s_eval;
        probitWorkingParam() = default;
        probitWorkingParam(arma::mat SS, arma::mat RR,arma::mat tXRR,arma::mat tXXx,arma::vec muu,int n){
            S = SS;
            tRR = SS*n;
            R = RR;
            tXR = tXRR;
            tXX = tXXx;
            mu = muu;
            s_eval = eig_sym(S);
        }
        void update(const arma::mat &Y,
                        const arma::mat &X,
                        const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        int n_rep, int nskp = 5);
};

// this will get us EM result of mu and the useful S matrix given Y and parameter at last iteration 
void probitWorkingParam::update(const arma::mat &Y,
                        const arma::mat &X,
                        const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        int n_rep, int nskp){

    int n = Y.n_rows;
    int q = Y.n_cols;
    arma::mat XB = X * B_t;
    XB.each_row() += mu_t.t();
    arma::mat C = arma::chol(Sigma_t);
    arma::vec Adiag(q);
    arma::vec b(q);
    arma::mat A(q,q);
    S.zeros(q,q);
    mu.zeros(q);
    R.zeros(n,q);
    tXR.zeros(q,q);
    s_eval.zeros(q);
    //arma::mat S(q,q,fill::zeros);
    //arma::vec mu(q,fill:zeros);

    for(int i = 0 ; i < n ; i++){
        Rcpp::checkUserInterrupt();
        Adiag = 2*arma::trans( Y.row(i) )-1;
        b = arma::trans(XB.row(i)) % Adiag;
        A = C.t();
        A.each_col() %= Adiag;
        Adiag = arma::solve(C.t(), Adiag-arma::trans(XB.row(i)));// this serves as the initial point
        LinearConstraints lincon(A,b,true);
        //Rcout << "      " << i <<"th sample sampling start" << endl;
        EllipticalSliceSampler sampler(n_rep + 1,lincon,nskp,Adiag);
        sampler.run();
        //Rcout << "      end" << endl;
        arma::mat resi = sampler.loop_state.samples;
        resi.shed_row(0);// we will remove the initial points
        resi = resi * C;
        R.row(i) = mean(resi);
        S += resi.t() * resi;
        resi.each_row() += mu_t.t();
        mu += arma::trans(arma::mean(resi)); 
    }

    S /= (n*n_rep);
    mu /= n;
    R.each_row() += mu_t.t();
    R.each_row() -= mu.t();
    tXR = X.t() * R;
    tRR = S*n;
    
    s_eval = eig_sym(S);

}


class cgprobitWorkingParam: public probitWorkingParam{
    public: 
        arma::mat S_Omega;
        arma::mat M;
        cgprobitWorkingParam() = default;
        cgprobitWorkingParam(arma::mat SS, 
            arma::mat RR,
            arma::mat tXRR,
            arma::mat tXXx,
            arma::vec muu,int n): probitWorkingParam(SS, RR, tXRR, tXXx, muu, n){
            S_Omega = SS;
            M = SS;
        }
        void update(const arma::mat &Y,
                        const arma::mat &X,
                        const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        const arma::mat &Omega_t,
                        int n_rep, int nskp = 5);

};

void cgprobitWorkingParam::update(const arma::mat &Y,
                        const arma::mat &X,
                        const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        const arma::mat &Omega_t,
                        int n_rep, int nskp){
    int n = Y.n_rows;
    int q = Y.n_cols;
    arma::mat XB = X * B_t;
    XB.each_row() += mu_t.t();
    M = XB.t() * XB / n;
    arma::mat C = arma::chol(Sigma_t);
    arma::vec Adiag(q);
    arma::vec b(q);
    arma::mat A(q,q);
    S.zeros(q,q);
    S_Omega.zeros(q,q);
    mu.zeros(q);
    R.zeros(n,q);
    tXR.zeros(q,q);
    s_eval.zeros(q);
    //arma::mat S(q,q,fill::zeros);
    //arma::vec mu(q,fill:zeros);
    arma::rowvec meani;

    for(int i = 0 ; i < n ; i++){
        meani = XB.row(i) * Sigma_t;
        Rcpp::checkUserInterrupt();
        Adiag = 2*arma::trans( Y.row(i) )-1;
        b = meani.t() % Adiag;
        A = C.t();
        A.each_col() %= Adiag;
        Adiag = arma::solve(C.t(), Adiag-meani.t());// this serves as the initial point
        LinearConstraints lincon(A,b,true);
        //Rcout << "      " << i <<"th sample sampling start" << endl;
        EllipticalSliceSampler sampler(n_rep + 1,lincon,nskp,Adiag);
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


void unitdiag(arma::mat & Sigma, arma::mat &Omega);

#endif