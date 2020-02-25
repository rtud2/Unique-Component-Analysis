#' rPCA: residual PCA
#'
#' `rPCA` returns the target data (scaled to the background) projected onto the orthogonal complement of the first `bg_components` principal components of the background
#' 
#' @param target Target dataset - dataset of interest (data.table)
#' @param bg Background dataset (data.table)
#' @param n_components number of Principal components to calculate for the target data, after being projected onto the orthogonal complement of the background
#' if NULL, then n_components chooses automatically based on finding the best linear spline with respect to squared-error.
#' @param bg_components number of background principal components used. Tuning parameter because this affects the span of the Orthogonal Complement.
#'  If vector, then list of rPC projections is returned for each item in the same order.
#'                   if NULL, then bg_components chooses automatically based on finding the best linear spline with respect to squared-error.
#' @param standardize (logical) default TRUE. Standardize the target to the colMean and sd of the background. Scale the background data.                   
#' @param return_all (logical) whether to return the background PCs and target projected on the Orthogonal Complement of the background
#' @param total_component (optional) argument to be passed into choose_pc if bg_components is unspecified.
#' @return  Data projected on the Orthogonal Complement contrastive principal components

rPCA = function(target, bg, n_components = NULL, bg_components = NULL, standardize = T, return_all = F, check = F, ...){
 
  
  if(check){
    #check that the dimensions of target and bg are the same
    if(ncol(target) != ncol(bg)){
      stop("Dimension mismatch: Please check target and bg have the same number of columns")
    }
    #check if matrices are all numeric
    if(target[, sum(!sapply(.SD, is.numeric))] > 0 | bg[, sum(!sapply(.SD, is.numeric))] > 0){
      stop("at least one column is not numeric")
    }
    
    #check if variance is zero in target and background.
    if(length(target[, {temp_sd = sapply(.SD, sd);
    which(temp_sd ==0 | !is.finite(temp_sd))}]) + 
    length(bg[, {temp_sd = sapply(.SD, sd);
    which(temp_sd ==0 | !is.finite(temp_sd))}])>0){
      stop("variance is zero in either target or background")
    }
      
  }
  
   if(standardize){
    target = scale(target, center = bg[, lapply(.SD, mean)], scale = bg[, lapply(.SD, sd)]);
    bg = bg[, lapply(.SD, scale)];
   }
  
  min_dim <- min(dim(bg))
  
  #convert to matrix types
  target <- data.matrix(target)
  bg <- data.matrix(bg)
  
  # Rotate the background
  if(length(bg_components) == 0){
    bg_svd_all <- irlba(bg)
    bg_components <- choose_pc(bg_svd_all$d)
    bg_component_vector <- bg_components
    warning(paste0("bg_components is NULL, auto-selecting ", bg_components, " components"))
    cat("\n");
  }else if(length(bg_components) > 0 & max(bg_components) < min_dim){
    bg_component_vector <- bg_components  
    bg_components <- max(bg_components)
  }else{
    stop("max(bg_components) exceeds the minimum dimension")
  }
  
  #calculating eigenvectors of the background data
  bg_svd <-irlba(bg, nv = bg_components)$v
  
  #projection onto the orthogonal complement
  returned_obj <- future_lapply(seq_along(bg_component_vector), function(zz){
    temp_bg_svd <- bg_svd[,1:bg_component_vector[zz]]
    oc_target <- target - tcrossprod(target %*% temp_bg_svd, temp_bg_svd)
    
    if(is.null(n_components)){
      oc_target_svd_all <- irlba(oc_target)
      n_components <- choose_pc(oc_target_svd_all$d)
      warning(paste0("n_components is NULL, auto-selecting ", n_components, " components"))
      cat("\n");
    }
    res_target_svd <- irlba(oc_target, nv = n_components)$v

    reduced_target <- oc_target %*% res_target_svd
    
    if(return_all){
      return(list("bg_components" = bg_component_vector[zz],"reduced_target" = reduced_target, "bg_svd" = bg_svd, "res_target_svd" = res_target_svd, "oc_target" = oc_target))
    }else{
      return(reduced_target)
    }    
    })
  return(returned_obj)
}

