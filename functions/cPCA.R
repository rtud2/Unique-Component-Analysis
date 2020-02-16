#' cPCA: Contrastive PCA
#'
#' Suppose you have two datasets, a `target` and a `background` where the `target` were a dataset containing information on cases
#'  and `background` contained information on controls, or uninteresting variation.  cPCA works directly on the covariance matrices
#'  and seeks to find the directions/rotation (eigenvector) that maximizes the variance explained in the `target` and minimize the
#'  variance explained in the `background`.
#'
#' @param target dataset - dataset of interest
#' @param bg Background dataset
#' @param n_components number of Principal components to calculate (default is 2)
#' @param alpha tuning parameter for how hard to weight minimizing the variance of the rotated background data
#' @param return_all default FALSE. (logical) whether top eigenvectors and the reduced background should be returned
#' @param standardize (logical) default TRUE. Standardize the target to the colMean and sd of the background. Scale the background data.
#' @return either return the reduced target data, or if `return_all=TRUE`, also return top eigenvectors and the reduced background


cPCA = function(target, bg, n_components = 2, alpha = 1, standardize = T, return_all = F){

  if(!is.matrix(target) | !is.matrix(bg)){
    target <- as.matrix(target)
    bg <- as.matrix(bg)
  }

  if(standardize){
    target = scale(target,center = colMeans(bg),scale = apply(bg, 2, sd));
    bg = scale(bg);
  }

  target_cov = cov(target);
  bg_cov = cov(bg);

  sigma = target_cov - alpha * bg_cov
  v_top <- eigen(sigma, symmetric = T)$vectors[,1:n_components]
  #reduced_target <- og_target %*% v_top
  reduced_target <- target %*% v_top
  reduced_bg <- bg %*% v_top

  if(return_all){
    return(list("reduced_target" = reduced_target, "reduced_bg" = reduced_bg, "vectors" = v_top))
  }else{
    return(list("reduced_target" = reduced_target))
  }
}

