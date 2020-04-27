#' @importFrom Rcpp evalCpp
#' @useDynLib uca, .registration = TRUE
NULL

#' multiple_score_calc
#' 
#' calculate the eigendecomposition via a product of matrices using SVD and QR decomposition. Much faster for large matrices.
#'  compute the lambda for the jth background 
#'
#' @param left a rectangular matrix with p columns
#' @param right a rectangular matrix with p rows
#' @param svd_right singular value decomposition of the right matrix
#' @param tau contrastive parameter
#' @return list of the largest eigenvalue, associated score, and tau
#' @importFrom methods as


multiple_score_calc <- function(left, right, svd_right, tau, B_focus){
  
  qr_left_U <- arma_qr(left %*% svd_right$u)
  RS_svd <- arma_svd( sweep(qr_left_U$R, 2, svd_right$d, FUN = "*") )
  u = qr_left_U$Q %*% RS_svd$u
  eigenvalues <- colSums(u * (left %*% (right %*% u))) #calculates diag(crossprod(u, left) %*% (right %*%u))
  top_eig_vals <- which.max(eigenvalues)
  
  return(
    list(score = as(arma_score(B_focus, matrix(u[,top_eig_vals], ncol = 1)), "numeric"),
         values = eigenvalues[top_eig_vals],
         tau = tau))
}


#' broken_svd
#' 
#' calculate SVD of a product of matrices by using svd and QR decompositions
#' 
#' @param left left side of a product
#' @param right right side of a product
#' @param nv number of unique components
#' @return top nv eigenvalues and associated eigenvectors
#' 
broken_svd = function(left, right, nv){
  svd_right <- arma_svd(right)
  qr_left_U <- arma_qr(left %*% svd_right$u)
  RS_svd <- arma_svd( sweep(qr_left_U$R, 2, svd_right$d, FUN = "*") )
  u = qr_left_U$Q %*% RS_svd$u
  
  eigenvalues <- colSums(u * (left %*% (right %*% u))) #calculates diag(crossprod(u, left) %*% (right %*%u))
  top_eig_vals <- order(eigenvalues, decreasing = T)[1:nv]
  list(values = eigenvalues[top_eig_vals],
       vectors = u[,top_eig_vals])
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
  
  right <- rbind(A , B )
  svd_right <- arma_svd(right)
  
  f_val = lapply(limit, function(z){
    left <- cbind(t(A ), - z * t(B ))
    multiple_score_calc(left, right, svd_right, z, B )})
  
  if(f_val[[1]]$score > 0){
    warning("Redundant Constraint: Lagrange Multiplier is negative. Setting lambda to 0 \n");
    return(f_val[[1]]);
  }else{
    
    for(iter in 1L:maxit){
      limit <- c(f_val[[1]]$tau, f_val[[2]]$tau)
      if( limit[2] - limit[1] < tol * limit[1]) break;
      if(iter == maxit) warning("maximum iteration reached: solution may not be optimal \n");
      
      tau = sum(limit)/2
      left <- cbind(t(A ), - tau * t(B ))
      tau_score = multiple_score_calc(left, right, svd_right, tau, B )
      
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
#' @importFrom future.apply future_sapply

magic_eigen_multiple = function(A, B, lambda, j, limit = c(0,20), maxit = 1E5, tol = 1E-6){
  
  #constants that don't really change if focused on j-th background
  right <- rbind(A, do.call(rbind, B[-j]), B[[j]])
  svd_right <- arma_svd(right)
  
  #checking bounds 
  f_val <- future_lapply(limit, function(zz) {
    left <- cbind( t(A), do.call(cbind, Map("*", -lambda[-j], lapply(B[-j], t))), -zz *t(B[[j]]));
    multiple_score_calc(left, right,svd_right, zz, B[[j]] ) })
  
  if(f_val[[1]]$score > 0){
    warning("Redundant Constraint: Lagrange Multiplier is negative. Setting lambda to 0 \n");
    return(f_val[[1]]);
  }else{
    for(iter in 1:maxit){
      limit <- c(f_val[[1]]$tau, f_val[[2]]$tau)
      if( limit[2] - limit[1] < tol * limit[1]) break;
      if(iter == maxit) warning("maximum iteration reached: solution may not be optimal \n");
      
      left <- cbind( t(A), do.call(cbind, Map("*", -lambda[-j], lapply(B[-j], t))), t(-sum(limit)/2 * B[[j]]))
      tau_score = multiple_score_calc(left, right, svd_right, sum(limit)/2, B[[j]] )
      
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
    for (j in seq_along(B )){
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
       valuess = bisection_j$values,
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
