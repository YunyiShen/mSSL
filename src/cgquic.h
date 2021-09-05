// cgquic.h
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef CGQUIC_H
#define CGQUIC_H

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <RcppArmadillo.h>

#define EPS (double(2.22E-16))

using namespace Rcpp;
using namespace arma;
using namespace std;

typedef struct
{
	int i;
	int j;
} ushort_pair_t;

// this updates the newton direction
void CoordinateDescentUpdate(
	int &q,						   // dimension
	const arma::mat &S,			   // Y^TY/n
	const arma::mat &Rho,		   // penalty
	const arma::mat &Omega,		   // precision matrix, X
	const arma::mat &Sigma,		   // covariance matrix, W=X^-1
	arma::mat &U,				   // DW
	const arma::mat &Q,			   // MW, new in gcquic
	arma::mat &D,				   // Newton directions, to be updated
	int i, int j,				   // looks like coordinates
	double &normD, double &diffD); // something related to the direction

// more or less from quic.h in mSSL, but make use of armadillo, pointers are difficult to read
cube cgquic(int &q,
			const arma::mat &S,	  // Y^TY/n
			const arma::mat &M,	  // u^Tu/n, new in gcquic
			const arma::mat &Rho, // hyper parameter about lasso
			double &tol,
			int &max_iter_quic);

#endif