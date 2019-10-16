rPCA_dz = function(target,
                   bg,
                   n_components = 2,
                   bg_components = 2 ## first k PCs of B explain variation in background
                   ) {
  X = data.matrix(target)
  #X = data.matrix(target)
  B = data.matrix(bg)
  v = svd(scale(B))$v # SVD on scaled background data, no difference here
  
  selected_v <- v[,-(1:bg_components)]
  newX = scale(X, center = colMeans(B), scale = apply(B, 2, sd)) %*%  selected_v # why scale to the background??
#  newX = X %*%  selected_v # why scale to the background??
  reduced_target = data.table(newX %*% svd(newX, nv = n_components)$v)
  return(reduced_target)  
}
