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


//[[Rcpp::export]]
arma::mat linconGauss_cpp(int n, 
                            arma::mat A, 
                            arma::vec b, 
                            arma::mat Sigma, 
                            arma::vec mu, 
                            arma::vec x_init,
                            int nskp=5){
    arma::mat C = chol(Sigma);
    b += A * mu;
    A = A * C.t();
    LinearConstraints lincon(A,b,true);
    x_init = arma::solve(C.t(), x_init-mu);
    
    EllipticalSliceSampler sampler(n,lincon,nskp,x_init);
    sampler.run();
    arma::mat res = sampler.loop_state.samples;
    res = res * C;
    res.each_row() += mu.t();

    return res;
}


class muSEM{
    public:
        arma::mat S;
        arma::vec mu;
        muSEM() = default;
        muSEM(arma::mat SS, arma::vec muu){
            S = SS;
            mu = muu;
        }
        void update(const arma::mat &Y,
                        const arma::mat &X,
                        const arma::vec &intercept,
                        const arma::mat &B,
                        const arma::mat &Sigma,
                        int n_rep, int nskp = 5);
};

// this will get us EM result of mu and the useful S matrix given Y and parameter at last iteration 
void muSEM::update(const arma::mat &Y,
                        const arma::mat &X,
                        const arma::vec &intercept,
                        const arma::mat &B,
                        const arma::mat &Sigma,
                        int n_rep, int nskp = 5){

    int n = Y.n_rows;
    int q = Y.n_cols;
    arma::mat XB = X * B;
    XB.each_row() += intercept.t();
    arma::mat C = arma::chol(Sigma);
    arma::vec Adiag(q);
    arma::vec b(q);
    arma::mat A(q,q);
    S.zeros(q,q);
    mu.zeros(q);
    //arma::mat S(q,q,fill::zeros);
    //arma::vec mu(q,fill:zeros);

    for(int i = 0 ; i < n ; i++){
        Adiag = 2*arma::trans( Y.row(i) )-1;
        b = arma::trans(XB.row(i)) % Adiag;
        A = C.t();
        A.each_col() %= Adiag;
        Adiag = arma::solve(C.t(), Adiag-arma::trans(XB.row(i)));// this serves as the initial point
        LinearConstraints lincon(A,b,true);
        EllipticalSliceSampler sampler(n_rep + 1,lincon,nskp,Adiag);
        sampler.run();
        arma::mat resi = sampler.loop_state.samples;
        resi.shed_row(0);// we will remove the initial points
        resi = resi * C;
        resi.each_row() += intercept.t();
        S += resi.t() * resi;
        mu += arma::trans(arma::mean(resi)); 
    }

    mu /= n;
    S /= (n*n_rep);
}
