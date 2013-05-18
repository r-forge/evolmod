kronecker.sum = function(matrix1, matrix2){
  
  if (!is.matrix(matrix1) || (nrow(matrix1) != ncol(matrix1))) 
    stop("\"matrix1\" must be a square matrix")
  
  if (!is.matrix(matrix2) || (nrow(matrix2) != ncol(matrix2))) 
    stop("\"matrix2\" must be a square matrix")
  
  dim1 = dim(matrix1)[1]
  dim2 = dim(matrix2)[2]
  
  ident.mat1 = diag(rep(1,dim1))
  ident.mat2 = diag(rep(1,dim2))
  
  ## %x% = Kronecker product
  
  return(ident.mat1%x%matrix2 + matrix1%x%ident.mat2)
}
