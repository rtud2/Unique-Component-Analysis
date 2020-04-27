// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenSandwich(const Eigen::Map<Eigen::MatrixXd>& v, Eigen::Map<Eigen::MatrixXd>& B){
  return Rcpp::wrap(v.transpose() * B * v);
}

