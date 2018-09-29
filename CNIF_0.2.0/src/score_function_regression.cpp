#include <RcppArmadillo.h>
//using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
double score_function_regression_with_fix(int node, const arma::uvec & parents, const arma::uvec & fixed, const arma::mat & data){
  arma::mat W = data.cols( join_vert(parents-1,fixed-1) );
  arma::mat X = data.cols( fixed );
  arma::mat Y = data.unsafe_col(node-1);
  arma::mat xtx = X.t() * X;
  arma::mat xtx_inv = pinv( (xtx.t() + xtx )/2.0 );
  arma::mat wtx = W.t() * X;
  arma::mat xtw = wtx.t();
  arma::mat wtxxtx_inv_xt =  wtx * xtx_inv * X.t();
  arma::mat tmp = wtxxtx_inv_xt * W;
  arma::mat delta = pinv( (tmp+tmp.t())/2 ) * wtxxtx_inv_xt * Y;
  
  arma::mat residual = pow(Y - W*delta,2);
  return( accu(residual) );
}

// [[Rcpp::export]]
double score_function_regression(int node, const arma::uvec & parents, const arma::mat & data){
  arma::mat X = data.cols( parents - 1 );
  arma::mat Y = data.unsafe_col(node - 1);
  arma::vec coef = solve(X,Y);  
  arma::mat residual = pow(Y - X*coef,2);
  return( accu(residual) );
}