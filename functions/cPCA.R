#' cPCA: Contrastive PCA
#'
#' Suppose you have two datasets, a `target` and a `background` where the `target` were a dataset containing information on cases
#'  and `background` contained information on controls, or uninteresting variation.  cPCA works directly on the covariance matrices
#'  and seeks to find the directions/rotation (eigenvector) that maximizes the variance explained in the `target` and minimize the
#'  variance explained in the `background`.
#'
#' @param target dataset - dataset of interest
#' @param bg Background dataset(s). Input multiple datasets as a list. 
#' @param n_components number of Principal components to calculate (default is 2)
#' @param bg_components number of principal components to keep in the background. (default is all) This induces a low rank approximation 
#' of the background covariance matrix
#' @param alpha tuning parameter for how hard to weight minimizing the variance of the rotated background data (default = 1)
#' if NULL, use Rayleigh Quotient and requires inversion of the `background` matrix. If positive semi-definite, a `fudge` is added.
#' @param return_all default FALSE. (logical) whether top eigenvectors and the reduced background should be returned
#' @param standardize (logical) default TRUE. Standardize the target to the colMean and sd of the background. Scale the background data.
#' @return either return the reduced target data, or if `return_all=TRUE`, also return top eigenvectors and the reduced background


cPCA = function(target, bg, n_components = 2, bg_components = ncol(bg), alpha = 1, standardize = F, return_all = F, ...){

  
  if(standardize){
    target = scale(target, center = bg[, lapply(.SD, mean)], scale = bg[, lapply(.SD, sd)]);
    bg = bg[, lapply(.SD, scale)];
  }else{
    target = scale(target, center = bg[,lapply(.SD, mean)], scale = F);
    bg = scale(bg, scale = F);
  }
  
  target <- data.matrix(target)
  bg <- data.matrix(bg)
  
  target_cov = crossprod(target)/nrow(target);
  
  #construct background covariance matrix
  if(between(bg_components, 1, ncol(bg), incbounds = F)){ # when bg_cov lies in low rank space
    top_bg <- irlba(bg, nv = bg_components)
    bg_cov <- tcrossprod(top_bg$v %*% diag(top_bg$d[1:bg_components]^2), top_bg$v)/nrow(bg)
  }else if(bg_components == ncol(bg)){
    bg_cov = crossprod(bg)/nrow(bg);
  }else{
    stop("bg_components is too large or negative")
  }
  
  if(length(alpha) == 0){
    #use rayleight quotient method, but make sure it's positive-definite
    fudge <- min(eigen(bg_cov, symmetric = T, only.values = T)$values)
    if(fudge <= 0){
      bg_cov <- bg_cov + diag(10*abs(fudge), nrow = nrow(bg_cov))
    }
    sigma = solve(bg_cov) %*% target_cov
  }else{ 
    sigma = target_cov - alpha * bg_cov
  }

  v_top <- eigen(sigma, symmetric = T)$vectors[,1:n_components]
  reduced_target <- target %*% v_top
  
  if(return_all){
      reduced_bg <- bg %*% v_top
      return(list("reduced_target" = reduced_target, "reduced_bg" = reduced_bg, "vectors" = v_top))
  }else{
    return(list("reduced_target" = reduced_target))
  }
}

