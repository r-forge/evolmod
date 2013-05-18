## generates a pairwise sequence alignment

sim.pair.align = function(prob.matrix, init.dist, aln.len){
  
  if (!is.matrix(prob.matrix) || (nrow(prob.matrix) != ncol(prob.matrix))) 
    stop("\"prob.matrix\" must be a square matrix")
  
  if (nrow(prob.matrix) != length(init.dist)) 
    stop("dimensions of prob matrix and init distribution do not match")
  
  state.space = c(1:length(init.dist))
  
  seq.table = matrix(rep(0,2*aln.len), nrow = 2, ncol = aln.len)
  
  seq.table[1,]= sample(state.space, aln.len, replace = T, init.dist)
  for (i in 1:aln.len){
    seq.table[2,i] = sample(state.space, 1, replace = T, prob.matrix[seq.table[1,i],])
  }
  
  return(seq.table)
}
