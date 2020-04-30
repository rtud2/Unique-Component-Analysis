#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
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

// [[Rcpp::export]]
SEXP broken_svd_cpp(const arma::mat & left, const arma::mat& right, const int & nv){
  arma::mat u, v; arma::vec d;
  arma::svd_econ(u,d,v,right, "left","dc");
  
  arma::mat q, r;
  arma::qr_econ(q, r, left * u);
  
  arma::mat rs_u, rs_v; arma::vec rs_d;
  
  //r.each_col() %= d; //multiply each column by corresponding element in d. no need to construct a diagonal matrix
  //arma::svd_econ(rs_u,rs_d,rs_v, r);
  arma::svd_econ(rs_u,rs_d,rs_v, r * diagmat(d));

  const arma::mat & evectors = q * rs_u;
  const arma::vec & evalues = arma::vectorise(arma::sum( evectors % (left*(right*evectors)), 0));
  
  const arma::uvec & idx = arma::sort_index(evalues,"descend");
  
  const arma::mat & out_evectors = evectors.cols(idx.subvec(0,nv-1));
  const arma::vec & out_evalue = evalues(idx.subvec(0, nv-1));
  List res; 
  res["values"] = wrap(out_evalue);
  res["vectors"] = wrap(out_evectors);
  return res;

  }

// [[Rcpp::export]]
SEXP multiple_score_calc_cpp(const arma::mat & left, const arma::mat& right, const arma::mat& right_u, const arma::vec & right_d){
  
  arma::mat q, r;
  arma::qr_econ(q, r, left * right_u);
  
  arma::mat rs_u, rs_v; arma::vec rs_d;
  
  //r.each_col() %= right_d; //multiply each column by corresponding element in d. no need to construct a diagonal matrix
  //arma::svd_econ(rs_u,rs_d,rs_v, r);
  arma::svd_econ(rs_u,rs_d,rs_v, r * diagmat(right_d));
  
  const arma::mat & evectors = q * rs_u;
  const arma::vec & evalues = arma::vectorise(arma::sum( evectors % (left*(right*evectors)), 0));
  
  const arma::uvec & idx = arma::sort_index(evalues,"descend");
  
  const arma::mat & out_evectors = evectors.cols(idx.head(1));
  const arma::vec & out_evalue = evalues(idx.head(1));

  List res; 
  res["vectors"] = out_evectors;
  res["values"] = wrap(out_evalue);
  return res;
  
}


