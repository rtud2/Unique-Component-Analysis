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
#' @importFrom RSpectra svds
#' 
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
    bg_svd_all <- svd(bg)
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
  bg_svd <-svds(bg, k = bg_components)$v
  
  #projection onto the orthogonal complement
  returned_obj <- lapply(seq_along(bg_component_vector), function(zz){
    temp_bg_svd <- bg_svd[,1:bg_component_vector[zz]]
    oc_target <- target - tcrossprod(target %*% temp_bg_svd, temp_bg_svd)
    
    if(is.null(n_components)){
      oc_target_svd_all <- svd(oc_target)
      n_components <- choose_pc(oc_target_svd_all$d)
      warning(paste0("n_components is NULL, auto-selecting ", n_components, " components"))
      cat("\n");
    }
    res_target_svd <- svds(oc_target, k = n_components)$v

    reduced_target <- oc_target %*% res_target_svd
    
    if(return_all){
      return(list("bg_components" = bg_component_vector[zz],"reduced_target" = reduced_target, "bg_svd" = bg_svd, "res_target_svd" = res_target_svd, "oc_target" = oc_target))
    }else{
      return(reduced_target)
    }    
    })
  return(returned_obj[[1]])
}

#' Choose PC Function
#' 
#' An automatic method of chooseing the number of principal components using the elbow method.
#' From the top 20 components, choose the best linear single knot spline. The linear spline with knot at `k` with the lowest MSE
#'  is the number of components we should choose. 
#' 
#' @param d the vector of diagonals from an SVD
#' @param total_component default is 20. The most number of components to include when fitting linear single knot splines 
#' @param return_all default FALSE.  If TRUE, then data.frame of mean squared error for the linear spline with knot for each principal component is outputted
#' @return either return the number optimal number of principal components or additionally, the data.frame of squared error for the linear spline 
#' with knot for each principal component is outputted
#' @importFrom splines bs

choose_pc <- function(d, total_component = 20, return_all = F){
  
  n_diag <- min(length(d), total_component) # capping at the number of components we look at to be the min of (n_components, total_component)
  upper <- n_diag - 1
  x <- 1:n_diag
  d <- d[1:n_diag]^2
  
  fit_list <- sapply(2:upper, function(xx) mean((d - lm(d ~ bs(x, degree = 1, knots = xx))$fitted.value)^2))
  res_dat <- data.frame(n_component = 2:upper, mse = fit_list)
  # best fitting linear spline at the elbow, means we take only the components before it
  n_comp <- res_dat[which.min(res_dat$mse), "n_component"] - 1
  if(return_all == T){
    return(list("n_component" = n_comp, "MSE" = res_dat))
  }else{
    return(c("n_component" = n_comp));
  }
  
}
