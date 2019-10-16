#' Choose PC Function
#' 
#' An automatic method of chooseing the number of principal components using the elbow method.
#' From the top 20 components, choose the best linear single knot spline. The linear spline with knot at `k` with the lowest MSE
#'  is the number of components we should choose. 
#' 
#' @param d the vector of diagonals from an SVD
#' @param return_all default FALSE.  If TRUE, then data.frame of squared error for the linear spline with knot for each principal component is outputted
#' @return either return the number optimal number of principal components or additionally, the data.frame of squared error for the linear spline 
#' with knot for each principal component is outputted

choose_pc <- function(d, return_all = F){
  
  n_diag <- min(length(d), 20) # capping at the number of components we look at to be the min of (n_components, 20)
  upper <- n_diag - 1
  x <- 1:n_diag
  d <- d[1:n_diag]^2
  
  fit_list <- lapply(2:upper, function(xx) lm(d ~ bs(x, degree = 1, knots = xx)))
  list_err <- lapply(fit_list, function(yy){
    plot_dat <- data.table(x, d, fit = predict(yy, newdata = data.frame(x = x)))
    plot_dat[, sum((d-fit)^2)]
  })
  res_dat <- data.table(n_component = 2:upper, sq_err = unlist(list_err))
  
  if(return_all == T){
    return(list("n_component" = res_dat[which.min(sq_err), n_component], "MSE" = res_dat))
  }else{
    return(c("n_component" = res_dat[which.min(sq_err), n_component]));
  }
  
}