#' rPCA: residual PCA
#'
#' `rPCA` returns the target data (scaled to the background) projected onto the orthogonal complement of the first `bg_components` principal components of the background
#' 
#' @param target Target dataset - dataset of interest
#' @param bg Background dataset
#' @param n_components number of Principal components to calculate for the target data, after being projected onto the orthogonal complement of the background
#' if NULL, then n_components chooses automatically based on finding the best linear spline with respect to squared-error.
#' @param bg_components number of background principal components used. Tuning parameter because this affects the span of the Orthogonal Complement
#'                   if NULL, then bg_components chooses automatically based on finding the best linear spline with respect to squared-error.
#' @param standardize (logical) default TRUE. Standardize the target to the colMean and sd of the background. Scale the background data.                   
#' @param return_all (logical) whether to return the background PCs and target projected on the Orthogonal Complement of the background
#' @param total_component (optional) argument to be passed into choose_pc if bg_components is unspecified.
#' @return  Data projected on the Orthogonal Complement contrastive principal components

rPCA = function(target, bg, n_components = NULL, bg_components = NULL, standardize = T, return_all = F, ...){
  if(!is.matrix(target) | !is.matrix(bg)){
    target <- as.matrix(target)
    bg <- as.matrix(bg)
  }
  og_target <- target
  
  if(standardize){
    target = scale(target,center = colMeans(bg),scale = apply(bg, 2, sd));
    bg = scale(bg);
  }
  # Rotate the background
  if(is.null(bg_components)){
    bg_svd_all <- svd(bg)
    bg_components <- choose_pc(bg_svd_all$d)
    warning(paste0("bg_components is NULL, auto-selecting ", bg_components, " components"))
    cat("\n");
    }
  bg_svd <-svd(bg, nv = bg_components)
  
  # Background Projection Matrix
  bg_projection <- tcrossprod(bg_svd$v)
  
  #projection onto the orthogonal complement
  oc_target <-  target %*% (diag(nrow = nrow(bg_projection)) - bg_projection)
  
  if(is.null(n_components)){
    oc_target_svd_all <- svd(oc_target)
    n_components <- choose_pc(oc_target_svd_all$d)
    warning(paste0("n_components is NULL, auto-selecting ", n_components, " components"))
    cat("\n");
  }
  res_target_svd <- svd(oc_target, nv = n_components)$v
  reduced_target <- oc_target %*% res_target_svd
  
  if(return_all){
    return(list("reduced_target" = reduced_target, "bg_svd" = bg_svd, "res_target_svd" = res_target_svd, "oc_target" = oc_target))
  }else{
    return(reduced_target)
  }
}

