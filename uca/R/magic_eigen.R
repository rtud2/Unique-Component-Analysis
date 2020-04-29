#' @importFrom Rcpp evalCpp
#' @useDynLib uca, .registration = TRUE
NULL


#' centeer
#' 
#' Convenient function to center the data, rather than typing `scale' or `sweep`
#' 
#' @param X a matrix
#' @return centered data matrix X
#' @export
#' 
center <- function(X){
  column_means <- colMeans(X)
  return(sweep(X, 2, column_means, "-"))
}


#' multiple_score_calc
#' 
#' calculate the eigendecomposition via a product of matrices using SVD and QR decomposition. Much faster for large matrices.
#'  compute the lambda for the jth background 
#'
#' @param left a rectangular matrix with p columns
#' @param right a rectangular matrix with p rows
#' @param svd_right a precalculated svd for the right matrix. saves time
#' @param tau contrastive parameter
#' @return list of the largest eigenvalue, associated score, and tau
#' @importFrom methods as


# multiple_score_calc <- function(left, right, tau, B_focus){
#   
#   tmp <- broken_svd_cpp(left, right , nv=1)
#   return(
#     list(score = as(arma_score(B_focus, tmp$vectors), "numeric"),
#          values =tmp$values,
#          tau = tau))
# }

multiple_score_calc <- function(left, right, svd_right, tau, B_focus){
  
 tmp_cpp <- multiple_score_calc_cpp(left, right, svd_right$u, svd_right$d)
  
  return(
    list(score = as(arma_score(B_focus, tmp_cpp$vectors), "numeric"),
         values = tmp_cpp$values,
         tau = tau))
}


#' bisection2
#' 
#' compute the UCA for single background using SVD and QR. good for bigger data where loading the covariance matrix is difficult
#' 
#' @param A a *Centered* n1xp data matrix
#' @param B a *Centered* n2xp data matrix
#' @param limit bounds for the lagrange multiplier.
#' @param maxit maximum iterations
#' @param tol tolerance for convergence criteria
#' @return list of tau, largest eigenvalue, and score
#' 
bisection2 = function(A, B, limit = c(0,20), maxit = 1E5L, tol = 1E-6){
  
  right <- rbind(A , B)
  svd_right <- arma_svd(right)
  t_A = t(A); 
  t_B = t(B);
  
  f_val <- vector(mode = "list", length = 2L)
  f_val[[1]] <- multiple_score_calc(left = t_A,
                                    right = A,
                                    svd_right = arma_svd(A),
                                    tau = 0,
                                    B_focus = B)
  f_val[[2]]$tau = 20
  
  if(f_val[[1]]$score > 0){
    warning("Redundant Constraint: Lagrange Multiplier is negative. Setting lambda to 0 \n");
    return(f_val[[1]]);
  }else{
    
    for(iter in 1L:maxit){
      limit <- c(f_val[[1]]$tau, f_val[[2]]$tau)
      if( limit[2] - limit[1] < tol * limit[1]) break;
      if(iter == maxit) warning("maximum iteration reached: solution may not be optimal \n");
      
      tau_score = multiple_score_calc(left = cbind(t_A, - (0.5*sum(limit)) * t_B),
                                      right = right,
                                      svd_right = svd_right,
                                      tau = (0.5*sum(limit)),
                                      B_focus = B )
      
      if(tau_score$score < 0){
        f_val[[1]] = tau_score
      }else{
        f_val[[2]] = tau_score
      }
    }  
    
    return(f_val[[ which.min(abs(c(f_val[[1]]$score, f_val[[2]]$score))) ]]) 
  }
}



#' magic_eigen_multiple
#' 
#' Solve for the optimal lagrange multiplier for the jth background. used only when multiple backgrounds exist.
#' 
#' @param A a *Centered* n1xp data matrix
#' @param B *Centered* a n2xp data matrix
#' @param lambda a j dimensional vector with possible lagrange multipliers for each background data matrix
#' @param j the specific background you're solving for
#' @param limit bounds for the lagrange multiplier.
#' @param maxit maximum iterations
#' @param tol tolerance for convergence criteria
#' @return list of tau, largest eigenvalue, and score

magic_eigen_multiple = function(A, B, lambda, j, limit = c(0,20), maxit = 1E5, tol = 1E-6){
  
  #constants that don't really change if focused on j-th background
  B_j <- B[[j]]
  t_B_j <- t( B_j )
  old_right <- rbind(A, do.call(rbind, B[-j]))
  right <- rbind(old_right, B_j)
  old_left <- cbind(t(A), do.call(cbind, Map("*", -lambda[-j], lapply(B[-j], t))))
  svd_right <- arma_svd(right)
  
  #checking bounds 
  f_val <- vector(mode = "list", length = 2L)
  f_val[[1]] <- multiple_score_calc(left = old_left,
                                    right = old_right,
                                    svd_right = arma_svd(old_right),
                                    tau = 0,
                                    B_focus = B_j)
  f_val[[2]]$tau = 20
  
  if(f_val[[1]]$score > 0){
    warning("Redundant Constraint: Lagrange Multiplier is negative. Setting lambda to 0 \n");
    return(f_val[[1]]);
  }else{
    for(iter in 1:maxit){
      limit <- c(f_val[[1]]$tau, f_val[[2]]$tau)
      if( limit[2] - limit[1] < tol * limit[1]) break;
      if(iter == maxit) warning("maximum iteration reached: solution may not be optimal \n");
      
      tau_score = multiple_score_calc(left = cbind( old_left, -0.5*sum(limit) * t_B_j),
                                      right = right,
                                      svd_right = svd_right,
                                      tau = 0.5*sum(limit),
                                      B_focus = B_j )
      
      if(tau_score$score < 0){
       f_val[[1]] <- tau_score
      }else{
       f_val[[2]] <- tau_score
      }
    }
    return(f_val[[ which.min(abs(c(f_val[[1]]$score, f_val[[2]]$score))) ]]) 
  }
}

#' bisection2.multiple
#' 
#' Solving for the unique components using SVD and QR in the instance of multiple backgrounds.
#' 
#' @param A *Centered* Target Data Matrix. n1 x p dimensions
#' @param B *Centered* List of k Background Data Matrix. n_k x p dimensions
#' @param lambda contrastive parameters if known. used to start algorithm. defaults to NULL
#' @param max_iter maximum number of iterations before giving up
#' @param tol convergence criteria
#' @return list of tau, largest eigenvalue, and score

bisection2.multiple <- function(A, B, lambda=NULL, max_iter=1E5L, tol = 1E-6, ...){
  
  #initialize starting point if one isn't supplied. greedy start
  if(length(lambda) == 0){
    #we use A and B here instead of divided b/c they divide in bisection2 function
    lambda = sapply(seq_along(B ), function(zz){bisection2(A, B[[zz]])$tau})
  }
  
  score = Inf; 
  
  #coordinate descent
  for(i in 1L:max_iter){
    old.score <- score
    for (j in seq_along(B)){
      #calculate the optimal lagrange multiplier for each background j
      bisection_j <- magic_eigen_multiple(A , B , lambda, j)
      lambda[j] <- bisection_j$tau
    }
    score <- sum(bisection_j$values, lambda)
    if( abs(old.score - score) < tol * abs(old.score)) break;
    #print(paste("iteration", i))
    }
  return(
    list(score = score,
       values = bisection_j$values,
       tau = lambda)
    )
}


 


# t_happy2 <- scale(t_happy, scale = F)
# t_neutral2 <- scale(t_neutral, scale = F)

###### for 100*136
# > microbenchmark::microbenchmark(magic_eigen(t_happy2, t_neutral2, 2, 5), eigs_sym(cov_happy - 2*cov_neutral, 5, "LA"), times = 10)
# Unit: seconds
# expr      min       lq     mean   median        uq       max neval cld
# magic_eigen(t_happy2, t_neutral2, 2, 5) 3.647757 3.667107 3.749784 3.696459  3.805642  4.017985    10  a 
# eigs_sym(cov_happy - 2 * cov_neutral, 5, "LA") 9.804725 9.869240 9.947251 9.920517 10.016461 10.162978    10   b


###### apparently cbinding transposes is easier than   rbinding then transpose.
# > microbenchmark::microbenchmark(t(rbind(A_divided, - lambda * B_divided)), cbind(t(A_divided), - lambda*t(B_divided)))
# Unit: milliseconds
# expr     min      lq     mean  median       uq      max neval cld
# t(rbind(A_divided, -lambda * B_divided)) 26.5603 26.8648 28.52919 27.0243 27.65665 142.5401   100   b
# cbind(t(A_divided), -lambda * t(B_divided)) 14.5453 16.4540 19.09077 16.5184 16.69970 140.5627   100  a 


#tmp_eigs <- magic_eigen(t_happy2, t_neutral2, 2, 5)
#tmp_eigs2 <- eigs_sym(cov_happy - 2*cov_neutral, 5, "LA")

#tmp_eigs$values - tmp_eigs2$values

#sign flips on eigenvectors
#sum(abs(tmp_eigs$vectors) - abs(tmp_eigs2$vectors))
