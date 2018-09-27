// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <mat.h>

using namespace Rcpp;

// shrinkage
mat shrinkage(const mat& a, double kappa);
RcppExport SEXP CNIF_shrinkage(SEXP aSEXP, SEXP kappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    __result = Rcpp::wrap(shrinkage(a, kappa));
    return __result;
END_RCPP
}
// covsel
mat covsel(const mat& d, const int max_iter, double alpha, double lambda, double rho);
RcppExport SEXP CNIF_covsel(SEXP dSEXP, SEXP max_iterSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    __result = Rcpp::wrap(covsel(d, max_iter, alpha, lambda, rho));
    return __result;
END_RCPP
}
// lasso
vec lasso(const mat& A, const vec& b, const int max_iter, double alpha, double lambda, double rho);
RcppExport SEXP CNIF_lasso(SEXP ASEXP, SEXP bSEXP, SEXP max_iterSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    __result = Rcpp::wrap(lasso(A, b, max_iter, alpha, lambda, rho));
    return __result;
END_RCPP
}
// ca_mean
arma::mat ca_mean(const arma::mat& x);
RcppExport SEXP CNIF_ca_mean(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    __result = Rcpp::wrap(ca_mean(x));
    return __result;
END_RCPP
}
// ca_var
arma::mat ca_var(const arma::mat& x);
RcppExport SEXP CNIF_ca_var(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    __result = Rcpp::wrap(ca_var(x));
    return __result;
END_RCPP
}
// center_scale
arma::mat center_scale(const arma::mat& x);
RcppExport SEXP CNIF_center_scale(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    __result = Rcpp::wrap(center_scale(x));
    return __result;
END_RCPP
}
// initial_sem_nofix
mat initial_sem_nofix(const mat& data, const uvec& variable, double rho, double lambda, double prec, int max_iter);
RcppExport SEXP CNIF_initial_sem_nofix(SEXP dataSEXP, SEXP variableSEXP, SEXP rhoSEXP, SEXP lambdaSEXP, SEXP precSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type variable(variableSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type prec(precSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    __result = Rcpp::wrap(initial_sem_nofix(data, variable, rho, lambda, prec, max_iter));
    return __result;
END_RCPP
}
// initial_sem
mat initial_sem(const mat& data, const uvec& fixed, const uvec& variable, double rho, double lambda, double prec, int max_iter);
RcppExport SEXP CNIF_initial_sem(SEXP dataSEXP, SEXP fixedSEXP, SEXP variableSEXP, SEXP rhoSEXP, SEXP lambdaSEXP, SEXP precSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type fixed(fixedSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type variable(variableSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type prec(precSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    __result = Rcpp::wrap(initial_sem(data, fixed, variable, rho, lambda, prec, max_iter));
    return __result;
END_RCPP
}
// score_function_regression_with_fix
double score_function_regression_with_fix(int node, const uvec& parents, const uvec& fixed, const mat& data);
RcppExport SEXP CNIF_score_function_regression_with_fix(SEXP nodeSEXP, SEXP parentsSEXP, SEXP fixedSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type parents(parentsSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type fixed(fixedSEXP);
    Rcpp::traits::input_parameter< const mat& >::type data(dataSEXP);
    __result = Rcpp::wrap(score_function_regression_with_fix(node, parents, fixed, data));
    return __result;
END_RCPP
}
// score_function_regression
double score_function_regression(int node, const uvec& parents, const mat& data);
RcppExport SEXP CNIF_score_function_regression(SEXP nodeSEXP, SEXP parentsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type parents(parentsSEXP);
    Rcpp::traits::input_parameter< const mat& >::type data(dataSEXP);
    __result = Rcpp::wrap(score_function_regression(node, parents, data));
    return __result;
END_RCPP
}
// get_cycles
map<int, vector<int> > get_cycles(const mat& adjMat);
RcppExport SEXP CNIF_get_cycles(SEXP adjMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type adjMat(adjMatSEXP);
    __result = Rcpp::wrap(get_cycles(adjMat));
    return __result;
END_RCPP
}
