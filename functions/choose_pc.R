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

choose_pc <- function(d, total_component = 20, return_all = F){
  
  n_diag <- min(length(d), total_component) # capping at the number of components we look at to be the min of (n_components, total_component)
  upper <- n_diag - 1
  x <- 1:n_diag
  d <- d[1:n_diag]^2
  
  fit_list <- sapply(2:upper, function(xx) mean((d - lm(d ~ bs(x, degree = 1, knots = xx))$fitted.value)^2))
  res_dat <- data.frame(n_component = 2:upper, mse = fit_list)
  
  if(return_all == T){
    return(list("n_component" = res_dat[which.min(res_dat$mse), "n_component"], "MSE" = res_dat))
  }else{
    return(c("n_component" = res_dat[which.min(res_dat$mse), "n_component"]));
  }
  
}
