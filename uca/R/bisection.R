#' score_calc
#' 
#' Calculate the derivative of Lagrangian.
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param tau lowerbound for root finding
#' @return list of tau (the lagrangian), eigenvector associated with tau, and derivative value of the lagrange multiplier
#' @importFrom methods as
#' @importFrom RSpectra eigs_sym
#' 
score_calc = function(A, B, tau){
  eigen_calc <- eigs_sym ( A - tau * B, 1L, "LA")
  return(list(#score = 1 - as(eigenSandwich(eigen_calc$vectors, B), "numeric"),
  score = 1 - as( crossprod(eigen_calc$vectors, B %*% eigen_calc$vectors), "numeric"),
              values = eigen_calc$values,
              tau = tau))
  }


#' bisection
#'
#' Use bisection method to find the optimal Lagrangian for solving unique component analysis (uca)
#' 
#' @param A Target Covariance Matrix
#' @param B Background Covariance Matrix. 
#' @param nv number of eigenvectors to use
#' @param limit upperbound of lagrange multiplier
#' @param maxit maxium number of iterations for the algorithm to run
#' @param tol tolerance for when to stop the algorithm
#' @return the final list of tau (the optimal lagrange multiplier), vector (eigenvector)
#'  associated with tau, value (eigenvalue), and score (derivative value of the lagrangian)

bisection = function(A, B, limit = 50, maxit = 1E5L, nv = 1, tol = 1E-6){
  
  f_val <- vector(mode = "list", length = 2L)
  f_val[[1]] <- score_calc(A, B, 0)
  og_upper_lim <- f_val[[2]]$tau <- limit

  if(f_val[[1]]$score >= 0){
    warning("Redundant Constraint: Lagrange Multiplier is negative. Setting lambda to 0 \n");
    return(f_val[[1]]);
  }else{
    
    for(iter in 1L:maxit){
      limit <- c(f_val[[1]]$tau, f_val[[2]]$tau)
      if( limit[2] - limit[1] < tol * limit[1]) break;
      if(iter == maxit) warning("maximum iteration reached: solution may not be optimal \n");
      
      tau_score = score_calc(A, B, sum(limit)/2)
      
      if(tau_score$score < 0){
        f_val[[1]] = tau_score
      }else{
        f_val[[2]] = tau_score
      }
    }  
    if(round(tau_score$tau) == og_upper_lim){
      warning(paste("Lagrange Multiplier is near upperbound. Consider increasing the upperbound. current limit is",og_upper_lim,"\n"))
    }
    return(f_val[[ which.min(abs(c(f_val[[1]]$score, f_val[[2]]$score))) ]]) 
  }
}


#' bisection.multiple
#' 
#' bisection method but for multiple background datasets. Use bisection method to find the optimal Lagrangian for solving UCA of each background dataset.
#' 
#' @param A Target Covariance Matrix
#' @param B list of background covariance(s).
#' @param lambda initial guess on what the Lagrange Multiplier should be. Default NULL
#' @param nv number of uca components to estimate
#' @param max_iter maximum iterations for coordinate descent, if tolerance is not reached. default 1E5
#' @param tol tolerance for stopping criteria of coordinate descent. default 1E-6
#' @param algo which algorithm to use. default algo == "bisection". If algo = "cd", L-BFGS-S optimization is used instead for coordinate descent. algo== "gd" than L-BFGS-S optimizes all lambda's simultaenously (gradient descent). 
#' @param ... other parameters to pass in to bisection(...)
#' @return for each background covariance matrix in B, return list of tau (the optimal lagrange multiplier), vectors(eigenvector)
#'  associated with largest value (eigenvalue)
#' @importFrom RSpectra eigs_sym

bisection.multiple = function(A, B, lambda=NULL, nv = 2L, max_iter = 1E5L, tol = 1E-6, algo = "bisection", ...){
  
  #initialize starting point if one isn't supplied
  if(length(lambda) == 0){
    lambda = c(0, sapply(seq_along(B)[-1], function(zz){bisection(A,B[[zz]])$tau}))
    }
  score = Inf; 
  if(algo == "bisection"){
    for(i in 1L:max_iter){
      old.score <- score
      for (j in seq_along(B)){
        bisection_j <- bisection(A - Reduce("+", Map("*", lambda[-j], B[-j])), B[[j]], ...) 
        lambda[j] = bisection_j$tau
      }
      score <- sum(bisection_j$values, lambda)
      if( abs(old.score - score) < tol * abs(old.score)) break;
    }
  }else if(algo == "cd"){
    for(i in 1L:max_iter){
      old.score <- score
      for (j in seq_along(B)){
        bisection_j <- optim_cd(A - Reduce("+", Map("*", lambda[-j], B[-j])), B[[j]], ...) 
        lambda[j] = bisection_j$tau
      }
      score <- sum(bisection_j$values, lambda)
      if( abs(old.score - score) < tol * abs(old.score)) break;
    }
  }else if(algo == "gd"){
        optim_gd_tmp <- optim_bfgs_gd(A = A, B = B, maxit = max_iter)
        lambda = optim_gd_tmp$tau
  }else{
    stop(paste("algo",algo," not recognized"))
  }
  
  
  dca <- eigs_sym(A - Reduce("+", Map("*", lambda, B)), nv, "LA")
 return(list(values = dca$values, vectors = dca$vectors, tau = lambda))
}
  

#' uca
#'
#' Run unique component analysis
#' 
#' @param A Target Data or Covariance Matrix
#' @param B list of background data or covariance matrices. 
#' @param nv number of uca components to estimate
#' @param method method used to calculate the uca values and vectors. Use method = 'data' when passing a n*p data matrix. Use method = 'cov' when passing in a covariance matrix 
#' @param center logical: default False. If False, data matrix A and B will not be centered
#' @param scale logical: default False. If True, will center and scale, regardless of what center variable is set to
#' @param algo algorithm to find lagrange multiplier. only takes values "bisection", "cd" (coordinate descent), and "gd" (gradient descent). For single background data, "cd" and "gd" are the same. default is "cd", but bisection exists for backwards compatibility. 
#' @param ... other parameters to pass in to bisection(...) and bisection.multiple(...).  (Default: limit=20, maxit=1E5L, max_iter = 1E5L, tol = 1E-6, algo = "bisection")
#' @return values(eigenvalues), vectors (eigenvectors), tau (Lagrange Multiplier) associated with unique component analysis
#' @importFrom RSpectra eigs_sym
#' @export

uca = function(A, B, nv = 2, method = "data", center = F, scale = F, algo = "cd", ...){
  if(sum(class(B) %in% c("list","matrix")) == 0 ){
    stop("B is not a list of matrix, matrices")
  }
  
  if(method == "cov"){
   # multiple background cov method 
    if(is.list(B) & length(B) > 1){
      if( sum((nrow(A) != nrow(A)), sapply(B , function(z){nrow(z) != ncol(z)})) > 0 ){
        stop("at least one input matrix is not square. make sure you've inputted a covariance matrix")
      }
      if(sum(sapply(lapply(B,dim), function(dims) all.equal(dims, dim(A)))) < length(B)){
        stop("at least one background dimension does not match target dimension")
      }
      
      bisection.multiple(A=A, B=B, nv=nv, algo=algo, ... )

      
    }else{
    # single background cov method
      if(is.list(B)) B = B[[1]]
      if((nrow(A) != nrow(B)) | nrow(A) != ncol(A) | nrow(B) != ncol(B)){
        stop("either A or B are not square, or don't have the same dimensions")
      }
      if(algo == "bisection"){
      tmp_res <- bisection(A=A, B=B, nv=nv, ...)
      }else if(algo == "cd" | algo == "gd"){
      tmp_res <- optim_cd(A=A, B=B, nv=nv, ...)  
      }else{
        stop(paste("algo", algo,"not regonized \n"))
      }
      res <- eigs_sym(A - tmp_res$tau*B, nv, "LA")
      return(list(values = res$values, vectors = res$vectors, tau = tmp_res$tau))
    }  
    
  }else if(method == "data"){
    
    if(is.list(B)){
      if(mean(sapply(B, ncol) == ncol(A)) < 1){
        stop("ncol(A) != ncol(B) in at least one element of list B")
      }
      if( sum(sapply(B, nrow) > ncol(A)) > 0){
        warning("Changing to method = 'cov' will possibly yield faster results.\n")
      }
    }else{
      # double checking dimensions
      nrows <- sapply(list(A,B), nrow)
      
      if(ncol(A) != ncol(B)){ #check the same number of variables
        stop("ncol(A) != ncol(B)")
      }
      if(sum(nrows > ncol(A)) > 0 ){ #if number of rows > columns, change method to cov 
        warning("Changing to method = 'cov' will possibly yield faster results.\n")
      }  
    }
  if( algo == "gd"){
      stop('algo == "gd" has not been implemented yet for method == "data"')
  }
    
    # scale data: single background
    if(scale == TRUE){
      A_divided = scale(A)/sqrt(nrow(A) - 1)
    }else if(center == TRUE){
      A_divided = center_f(A)/sqrt(nrow(A) - 1)
    }else{
      A_divided = A/sqrt(nrow(A) - 1)      
    }
    # scale data: Multiple backgrounds
    if(is.list(B) & length(B) > 1){
      #run multi-background
      if(scale == TRUE){
        B_divided <- Map(function(z){scale(z)/sqrt(nrow(z) - 1)}, B)
      }else if(center == TRUE){
        B_divided <- Map(function(z){center_f(z)/sqrt(nrow(z) - 1)}, B)
      }else{
        B_divided <- Map(function(z){z/sqrt(nrow(z) - 1)}, B)  
      }
      
      bisection2.multiple(A=A_divided, B=B_divided, nv = nv, algo = algo, ... )
      
      }else{
        #run single background
      if(is.list(B)) B = B[[1]]
      if(scale == TRUE){
        B_divided = scale(B)/sqrt(nrow(B) - 1)
      }else if(center == TRUE){
        B_divided = center_f(B)/sqrt(nrow(B) - 1)
      }else{
        B_divided = B/sqrt(nrow(B) - 1)    
      }
      if(algo == "bisection"){
        tmp_res <- bisection2(A=A_divided, B=B_divided, ...)  
      } else if(algo == "cd"){
        tmp_res <- optim_cd2(A = A_divided, B = B_divided, ...)
      }else{
        stop(paste("algo", algo,"not regonized \n"))
      }
      
      
      #calculate the svd
      left <- cbind(t(A_divided), - tmp_res$tau * t(B_divided))
      right <- rbind(A_divided, B_divided)
     
      final_res <- broken_svd_cpp(left, right, nv)
      return(list(values = final_res$values, vectors = final_res$vectors, tau = tmp_res$tau))
    }  
  }else{
    stop("Method is not 'cov' or 'data'")
  }
  
}
