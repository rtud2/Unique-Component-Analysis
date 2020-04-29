## minimize lambda_max(U - l * V) + l over scalar li
opt_l = function(U, V, max_l, min_l, maxit, tol) {
  l = (max_l + min_l) / 2
  for (it in 1:maxit) {
    if ((max_l - min_l) / min_l <= tol) {
      break
    } else {
      pca = partial_eigen(U - l * V, n = 1)
      v = pca$vectors
      ## derivative
      d = 1 - crossprod(v, V) %*% v
      ## update guess and bounds
      if (d >= 0) {
        max_l = l
        l = (min_l + l) / 2
      } else {
        min_l = l
        l = (l + max_l) / 2
      }
    }
  }
  maxit_flag = FALSE
  ## was maxit reached?
  if (it == maxit) {
    maxit_flag = TRUE
  }
  ## L = lagrangian, used later for checking convergence
  return(list(l = l, L = pca$values, maxit_flag = maxit_flag))
}

## k = number of eigenvectors
dcomp = function(tg, bg, k,
                 max_l = 100, min_l = 0,
                 maxit = 100, tol = 1e-6) {
  ## initial values
  l = rep((max_l + min_l) / 2, length(bg))
  L_old = 0
  L_new = Inf
  ## coordinate descent
  for (it in 1:maxit) {
    if (abs(L_old - L_new) < tol * L_old) {
      break
    } else {
      ## update lagrangian
      L_old = L_new
      ## loop over coordinates
      for (i in 1:length(bg)) {
        ## calculate U and V matrices
        U = tg
        for (j in setdiff(1:length(bg), i)) {
          U = U - l[j] * bg[[j]]
        }
        V = bg[[i]]
        ## minimize l[i]
        opt_l_i = opt_l(U, V, max_l, min_l, maxit, tol)
        l[i] = opt_l_i$l
        ## did this optimization reach maximum iterations?
        maxit_flag = opt_l_i$maxit_flag
      }
      ## new lagrangian is the lagrangian of the last opt_l
      ## run plus the sum of the lambdas
      L_new = opt_l_i$L + sum(l)
    }
  }
  print(l)
  if (it == maxit || maxit_flag) {
    warning("maxit reached")
  }
  ## calculate both eigenvectors using last U and V
  pe = partial_eigen(U - l[i] * V, n = k)
  v = pe$vectors
  return(v)
}
