// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mpSSL_dpe
List mpSSL_dpe(arma::mat X, arma::mat Y, List lambdas, List xis, arma::vec theta_hyper_params, arma::vec eta_hyper_params, int diag_penalty, int max_iter, double eps, int s_max_condition, int obj_counter_max, int verbose, int n_rep, int nskp);
RcppExport SEXP _mpSSL_mpSSL_dpe(SEXP XSEXP, SEXP YSEXP, SEXP lambdasSEXP, SEXP xisSEXP, SEXP theta_hyper_paramsSEXP, SEXP eta_hyper_paramsSEXP, SEXP diag_penaltySEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP s_max_conditionSEXP, SEXP obj_counter_maxSEXP, SEXP verboseSEXP, SEXP n_repSEXP, SEXP nskpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< List >::type xis(xisSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_hyper_params(theta_hyper_paramsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_hyper_params(eta_hyper_paramsSEXP);
    Rcpp::traits::input_parameter< int >::type diag_penalty(diag_penaltySEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type s_max_condition(s_max_conditionSEXP);
    Rcpp::traits::input_parameter< int >::type obj_counter_max(obj_counter_maxSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< int >::type nskp(nskpSEXP);
    rcpp_result_gen = Rcpp::wrap(mpSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params, diag_penalty, max_iter, eps, s_max_condition, obj_counter_max, verbose, n_rep, nskp));
    return rcpp_result_gen;
END_RCPP
}
// mpSSL_dcpe
List mpSSL_dcpe(arma::mat X, arma::mat Y, List lambdas, List xis, arma::vec theta_hyper_params, arma::vec eta_hyper_params, int diag_penalty, int max_iter, double eps, int verbose, int n_rep, int nskp);
RcppExport SEXP _mpSSL_mpSSL_dcpe(SEXP XSEXP, SEXP YSEXP, SEXP lambdasSEXP, SEXP xisSEXP, SEXP theta_hyper_paramsSEXP, SEXP eta_hyper_paramsSEXP, SEXP diag_penaltySEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP n_repSEXP, SEXP nskpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< List >::type xis(xisSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_hyper_params(theta_hyper_paramsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_hyper_params(eta_hyper_paramsSEXP);
    Rcpp::traits::input_parameter< int >::type diag_penalty(diag_penaltySEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< int >::type nskp(nskpSEXP);
    rcpp_result_gen = Rcpp::wrap(mpSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params, diag_penalty, max_iter, eps, verbose, n_rep, nskp));
    return rcpp_result_gen;
END_RCPP
}
// cgpSSL_dpe
List cgpSSL_dpe(arma::mat X, arma::mat Y, List lambdas, List xis, arma::vec theta_hyper_params, arma::vec eta_hyper_params, int diag_penalty, int max_iter, double eps, int s_max_condition, int obj_counter_max, int verbose, int n_rep, int nskp);
RcppExport SEXP _mpSSL_cgpSSL_dpe(SEXP XSEXP, SEXP YSEXP, SEXP lambdasSEXP, SEXP xisSEXP, SEXP theta_hyper_paramsSEXP, SEXP eta_hyper_paramsSEXP, SEXP diag_penaltySEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP s_max_conditionSEXP, SEXP obj_counter_maxSEXP, SEXP verboseSEXP, SEXP n_repSEXP, SEXP nskpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< List >::type xis(xisSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_hyper_params(theta_hyper_paramsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_hyper_params(eta_hyper_paramsSEXP);
    Rcpp::traits::input_parameter< int >::type diag_penalty(diag_penaltySEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type s_max_condition(s_max_conditionSEXP);
    Rcpp::traits::input_parameter< int >::type obj_counter_max(obj_counter_maxSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< int >::type nskp(nskpSEXP);
    rcpp_result_gen = Rcpp::wrap(cgpSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params, diag_penalty, max_iter, eps, s_max_condition, obj_counter_max, verbose, n_rep, nskp));
    return rcpp_result_gen;
END_RCPP
}
// unitdiag
void unitdiag(arma::mat& Sigma, arma::mat& Omega);
RcppExport SEXP _mpSSL_unitdiag(SEXP SigmaSEXP, SEXP OmegaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Omega(OmegaSEXP);
    unitdiag(Sigma, Omega);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mpSSL_mpSSL_dpe", (DL_FUNC) &_mpSSL_mpSSL_dpe, 14},
    {"_mpSSL_mpSSL_dcpe", (DL_FUNC) &_mpSSL_mpSSL_dcpe, 12},
    {"_mpSSL_cgpSSL_dpe", (DL_FUNC) &_mpSSL_cgpSSL_dpe, 14},
    {"_mpSSL_unitdiag", (DL_FUNC) &_mpSSL_unitdiag, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_mpSSL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
