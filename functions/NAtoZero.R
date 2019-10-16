## convert NA values to zero
NAtoZero = function(mat){
  mat[is.na(mat)] <- 0
  return(mat)
}