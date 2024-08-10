#ifndef WORKINGPARAM_PROBITPARAM_H
#define WORKINGPARAM_PROBITPARAM_H


// [[Rcpp::depends(RcppArmadillo)]]
#include <linconGaussR.h> 
#include "wparam.h"

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace linconGaussR;

namespace workingparam{

class probitWorkingParam : public WorkingParam{
    public:
        probitWorkingParam() = default;
        probitWorkingParam(arma::mat X, arma::mat Y): WorkingParam(X,Y) {}
        inline void update(const arma::vec &mu_t,
                        const arma::mat &B_t,
                        const arma::mat &Sigma_t,
                        int n_rep, int nskp = 5);
        inline void postprocessing(
                        arma::mat &B, 
                        arma::mat & Sigma, arma::mat &Omega){
            unitdiag(Sigma,Omega);
        }
};

// this will get us EM result of mu and the useful S matrix given Y and parameter at last iteration 
inline void probitWorkingParam::update(
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



class probitcontWorkingParam : public WorkingParam{
  public:
    int binidxend; // last index for binary responses
    probitcontWorkingParam() = default;
    probitcontWorkingParam(arma::mat X, arma::mat Y, int binidxend): WorkingParam(X,Y), binidxend(binidxend) {}
    inline void update(const arma::vec &mu_t,
                     const arma::mat &B_t,
                     const arma::mat &Sigma_t,
                     int n_rep, int nskp = 5);
    inline void postprocessing(
        arma::mat &B, 
        arma::mat & Sigma, 
        arma::mat &Omega){
        subunitdiag(Sigma, Omega, binidxend);
    }
};

// this will get us EM result of mu and the useful S matrix given Y and parameter at last iteration 
inline void probitcontWorkingParam::update(
    const arma::vec &mu_t,
    const arma::mat &B_t,
    const arma::mat &Sigma_t,
    int n_rep, int nskp){
  
    int q_tot = Y.n_cols; // total number of things
    int q = binidxend + 1; // number of binary things
    int n = Y.n_rows;
    
    arma::mat XB = X * B_t;
    XB.each_row() += mu_t.t();
    
    // usual continuous updates
    R.tail_cols(q_tot - q) = Y.tail_cols(q_tot - q)-XB.tail_cols(q_tot - q);
  
    
    
    // binary stuff
    XB = XB.cols(0, binidxend); 
    
    // to construct the linear constraints by transforming constraints with Chol
    arma::mat C = arma::chol(Sigma_t.submat(0,0,binidxend,binidxend)); // take only submatrix of Sigma
    arma::vec Adiag(q);
    arma::vec b(q);
    arma::mat A(q,q);
    // some temp matrices
    arma::mat Stmp, Rtmp, tXRtmp;
    arma::vec mutmp;
    Stmp.zeros(q,q);
    mutmp.zeros(q);
    Rtmp.zeros(n,q);
    s_eval.zeros(q_tot);
    
  
    for(int i = 0 ; i < n ; i++){
      Rcpp::checkUserInterrupt();
      Adiag = 2*arma::trans( Y.row(i) )-1; // this makes \pm 1
      b = arma::trans(XB.row(i)); 
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
      Rtmp.row(i) = mean(resi);
      Stmp += resi.t() * resi;
      resi.each_row() += mu_t(arma::span(0, binidxend)).t();
      mutmp += arma::trans(arma::mean(resi)); 
    }
  
    Stmp /= (n*n_rep);
    mutmp /= n;
    
    S.submat(0,0,binidxend,binidxend) = Stmp;
    R.submat(0,0,binidxend,binidxend) = Rtmp;
    mu(arma::span(0, binidxend)) = mutmp;
    
    // this is a bit wasteful since only part of mu is updated 
    R.each_row() += mu_t.t();
    R.each_row() -= mu.t();
    tXR = X.t() * R;
    tRR = S*n;
  
    s_eval = eig_sym(S);
  
}



}
#endif