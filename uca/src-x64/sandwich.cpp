#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
SEXP arma_qr(const arma::mat & X){
  arma::mat q, r;
  arma::qr_econ(q, r, X);
  List res; res["Q"] = wrap(q); res["R"] = wrap(r);
  return res;
}

// [[Rcpp::export]]
SEXP arma_svd(const arma::mat & X){
  arma::mat u, v; arma::vec d;
  arma::svd_econ(u,d,v,X,"left","dc");
  List res; res["u"] = wrap(u); res["d"] = wrap(d);
  return res;
}

// [[Rcpp::export]]
SEXP arma_score(const arma::mat & A, const arma::mat & B){
  const arma::mat & X = A * B;
  return wrap(1.0 - X.t() * X);
}