#' cPCA: Contrastive PCA
#'
#' Suppose you have two datasets, a `target` and a `background` where the `target` were a dataset containing information on cases
#'  and `background` contained information on controls, or uninteresting variation.  cPCA works directly on the covariance matrices
#'  and seeks to find the directions/rotation (eigenvector) that maximizes the variance explained in the `target` and minimize the
#'  variance explained in the `background`.
#'
#' @param A covariance of interest
#' @param bg Background dataset(s). Input multiple datasets as a list. 
#' @param nv number of Principal components to calculate (default is 2)
#' @param bg_components number of principal components to keep in the background. (default is all) This induces a low rank approximation 
#' of the background covariance matrix
#' @param alpha tuning parameter for how hard to weight minimizing the variance of the rotated background data (default = 1)
#' if NULL, use Rayleigh Quotient and requires inversion of the `background` matrix. If positive semi-definite, a `fudge` is added.
#' @param return_all default FALSE. (logical) whether top eigenvectors and the reduced background should be returned
#' @param standardize (logical) default TRUE. Standardize the target to the colMean and sd of the background. Scale the background data.
#' @return either return the reduced target data, or if `return_all=TRUE`, also return top eigenvectors and the reduced background
#' @importFrom RSpectra eigs_sym
#' @export
#' 

#cPCA = function(targ, bg, nv = 2, bg_components = ncol(bg), alpha = NULL, ...){
  
cPCA = function(A, B, alpha = NULL, nv = 2, ...){
  
  #construct background covariance matrix
  # if(between(bg_components, 1, ncol(bg), incbounds = F)){ # when bg_cov lies in low rank space
  #   top_bg <- irlba(bg, nv = bg_components)
  #   bg_cov <- tcrossprod(top_bg$v %*% diag(top_bg$d[1:bg_components]^2), top_bg$v)/nrow(bg)
  # }else if(bg_components == 1){
  #   top_bg <- irlba(bg, nv = bg_components)
  #   bg_cov <- tcrossprod(top_bg$v)*top_bg$d[1]^2/nrow(bg)
  # }else if(bg_components == ncol(bg)){
  #   bg_cov = cov(bg);
  # }else{
  #   stop("bg_components is too large or negative")
  # }
  
  if(length(alpha) == 0){
    #use rayleigh quotient method, but make sure it's positive-definite
    fudge <- eigs_sym(B, 1, which = "SA")$values
    if(fudge <= 0){
      diag(B) <- diag(B) + 10*abs(fudge)
    }
    sigma = solve(B) %*% A
  }else{ 
    sigma = A - alpha * B
  }

  res <- eigs_sym(sigma, nv, "LA")
  
  return(list(values = res$values, vectors = res$vectors, alpha = alpha))
}

