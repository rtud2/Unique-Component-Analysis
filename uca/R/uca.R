#' @aliases uca-package
"_PACKAGE"

#' Unique Component Analysis
#'
#' Run unique component analysis
#' @param A Target Data or Covariance Matrix
#' @param B list of background data or covariance matrices.
#' @param nv number of uca components to estimate
#' @param method method used to calculate the uca values and vectors.
#' Use method = 'data' when passing a n x p data matrix. Recommended when p > n
#' Use method = 'cov' when passing in a covariance matrix.
#' @param center logical: default FALSE. If TRUE, data matrix A & B are centered
#' @param scale logical: default FALSE. If TRUE, data matrix A & B are centered
#' and scaled, regardless of whether `center == TRUE` or `center == FALSE`
#' @param algo algorithm to find lagrange multiplier(s). valid algorithms are
#' "bisection", "cd" (coordinate descent), and "gd" (gradient descent). For
#' single background data, "cd" and "gd" are the same. default is "cd", but
#' bisection exists for backwards compatibility.
#' @param ... additional parameters for `bisection_cov()`, `optim_cov_cd()`
#' when `algo = "bisection"`:
#' \itemize{
#'  \item `limit` upper bound of lagrange multiplier
#'  \item `maxit` maximum number of iterations for algorithm to run
#'  \item `tol` tolerance to stop the algorithm
#'  \item `max_iter` maximum iteration for coordinate descent (multi-background)
#' when `algo = "cd"` or `algo = "gd"`:
#'  * `maxit` maximum number of iterations for bfgs algorithm to run
#' @return list of three elements:
#' \itemize{
#'  \item values (eigenvalues)
#'  \item vectors (eigenvectors) of unique component analysis
#'  \item tau optimal Lagrange Multiplier(s) associated with uca
#'  }
#' @importFrom RSpectra eigs_sym
#' @export

uca <- function(A, B, nv = 2, method = "data", center = FALSE, scale = FALSE,
                algo = "cd", ...) {
  if (sum(class(B) %in% c("list", "matrix")) == 0) {
    stop("B is not a list of matrix, matrices")
  }
  if (method == "cov") {
   # multiple background cov method
    if (is.list(B) & length(B) > 1) {
      if (sum((nrow(A) != nrow(A)),
              sapply(B, function(z) {
                         nrow(z) != ncol(z)
            })) > 0) {
        stop("at least one input matrix is not square.
             make sure you've inputted a covariance matrix")
      }

      if (sum(sapply(lapply(B, dim), function(dims) {
                         all.equal(dims, dim(A))
            })) < length(B)) {
        stop("at least one background dimension 
             does not match target dimension")
      }

      cov_multiple(A = A, B = B, nv = nv, algo = algo, ...)
    } else {
    # single background cov method
      if (is.list(B)) {
          B <- B[[1]]
      }
      if ((nrow(A) != nrow(B)) | nrow(A) != ncol(A) | nrow(B) != ncol(B)) {
        stop("either A or B are not square,
             or don't have the same dimensions")
      }
      if (algo == "bisection") {
        tmp_res <- bisection_cov(A = A, B = B, nv = nv, ...)
      } else if (algo == "cd" | algo == "gd") {
        tmp_res <- optim_cov_cd(A = A, B = B, nv = nv, ...)
      } else {
        stop(paste("algo", algo, "not regonized \n"))
      }
      res <- eigs_sym(A - tmp_res$tau * B, nv, "LA")
      return(list(values = res$values,
                  vectors = res$vectors,
                  tau = tmp_res$tau))
    }
  } else if (method == "data") {
    if (is.list(B)) {
      if (mean(sapply(B, ncol) == ncol(A)) < 1) {
        stop("ncol(A) != ncol(B) in at least one element of list B")
      }
      if (sum(sapply(B, nrow) > ncol(A)) > 0) {
        warning("Changing to method = 'cov'
                will possibly yield faster results.\n")
      }
    } else {
      # double checking dimensions
      nrows <- sapply(list(A, B), nrow)

      if (ncol(A) != ncol(B)) {
      # check the same number of variables
        stop("ncol(A) != ncol(B)")
      }
      if (sum(nrows > ncol(A)) > 0) {
      # if number of rows > columns, change method to cov
        warning("Changing to method = 'cov'
                will possibly yield faster results.\n")
      }
    }
  if (algo == "gd") {
      stop('algo == "gd" has not been implemented yet
           for method == "data"')
  }

    # scale data: single background
    if (scale == TRUE) {
      A_divided <- scale(A) / sqrt(nrow(A) - 1)
    } else if (center == TRUE) {
      A_divided <- center_f(A) / sqrt(nrow(A) - 1)
    } else {
      A_divided <- A / sqrt(nrow(A) - 1)
    }
    # scale data: Multiple backgrounds
    if (is.list(B) & length(B) > 1) {
      #run multi-background
      if (scale == TRUE) {
        B_divided <- Map(function(z) {
                             scale(z) / sqrt(nrow(z) - 1)},
                             B)
      } else if (center == TRUE) {
        B_divided <- Map(function(z) {
                             center_f(z) / sqrt(nrow(z) - 1)},
                             B)
      } else {
        B_divided <- Map(function(z) {
                             z / sqrt(nrow(z) - 1)},
                             B)
      }

      data_multiple(A = A_divided,
                    B = B_divided,
                    nv = nv,
                    algo = algo,
                    ...)

      } else {
      #run single background
      if (is.list(B)) {
        B <- B[[1]]
      }
      if (scale == TRUE) {
        B_divided <- scale(B) / sqrt(nrow(B) - 1)
      } else if (center == TRUE) {
        B_divided <- center_f(B) / sqrt(nrow(B) - 1)
      } else {
        B_divided <- B / sqrt(nrow(B) - 1)
      }
      if (algo == "bisection") {
        tmp_res <- bisection_data(A = A_divided,
                                  B = B_divided,
                                  ...)
      } else if (algo == "cd") {
        tmp_res <- optim_data_cd(A = A_divided,
                                 B = B_divided,
                                 ...)
      } else {
        stop(paste("algo", algo, "not regonized \n"))
      }
      #calculate the svd
      left <- cbind(t(A_divided), - tmp_res$tau * t(B_divided))
      right <- rbind(A_divided, B_divided)

      final_res <- broken_svd_cpp(left, right, nv)
      return(list(values = final_res$values,
                  vectors = final_res$vectors,
                  tau = tmp_res$tau))
    }
  } else {
    stop("Method is not 'cov' or 'data'")
  }
}
