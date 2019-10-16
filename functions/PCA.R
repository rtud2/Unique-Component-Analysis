#' PCA: Principal Component Analysis
#' 
#' Suppose you have two datasets, a `target` and a `background` where the `target` were a dataset containing information on cases
#'  and `background` contained information on controls, or uninteresting variation.  cPCA works directly on the covariance matrices 
#'  and seeks to find the directions/rotation (eigenvector) that maximizes the variance explained in the `target` and minimize the 
#'  variance explained in the `background`.
#' 
#' @param target dataset - dataset of interest
#' @param n_components number of Principal components to calculate (default is 2)
#' @param alpha tuning parameter for how hard to weight minimizing the variance of the rotated background data
#' @param standardize (logical) default TRUE. Standardize the target
#' @param return_all default FALSE. (logical) whether top eigenvectors should be returned
#' @return either data projected on the principal components, or if `return_all=TRUE`, also return top eigenvectors


PCA = function(target, n_components = 2, standardize = T, return_all = F){
  if(!is.matrix(target)){
    target <- as.matrix(target)
  }
  if(standardize){
    target = scale(target);
  }
  target_cov = cov(target);
  v_top <- svd(target_cov, nv = n_components)$v
  reduced_target <- target %*% v_top 
  if(return_all){
    return(list("reduced_target"=reduced_target, "vectors"=v_top))
  }else{
    return(list("reduced_target"=reduced_target))  
  }
}