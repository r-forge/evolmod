as.eigen = function (x, ...) UseMethod("as.eigen")

as.eigen.default = function(x){
  markov.eigen = eigen(x)
  inv.eigen = solve(markov.eigen$vectors)
  
  markov.chain = list()
  markov.chain[["rate.matrix"]] = x
  markov.chain[["values"]] = markov.eigen$values
  markov.chain[["vectors"]] = markov.eigen$vectors
  markov.chain[["invvectors"]] = inv.eigen
  class(markov.chain) = "eigen"
  return(markov.chain)
}

as.eigen.revmc = function(x){
#  if (!("revmc" %in% class(x)))
#    stop("Error: object \"revmc.object\" is not of class \"revmc\"")
  
  diag.stat = diag(x$stationary^0.5)
  diag.stat.inv = diag(x$stationary^-0.5)
  
  aux.matrix = diag.stat%*%x$rate.matrix%*%diag.stat.inv
  
  aux.eigen = eigen(aux.matrix, symmetric = T)
  
  revmc.eigen = list(x$rate.matrix,
                     x$stationary, aux.eigen$values,
                     diag.stat.inv%*%aux.eigen$vectors,
                     t(aux.eigen$vectors)%*%diag.stat)
  names(revmc.eigen) = c("rate.matrix", "stationary",
                         "values", "vectors", "invvectors")
  
  class(revmc.eigen) = c("revmc", "eigen")
  
  return(revmc.eigen)  
}

matexp = function(rate.eigen, interval.len){
  if (!("revmc" %in% class(rate.eigen)))
    stop("Error: object \"rate.eigen\" is not of class \"revmc\"")
  
  if (!("eigen" %in% class(rate.eigen)))
    stop("Error: object \"rate.eigen\" is not of class \"eigen\"")
  
  prob.matrix = NULL
  
  ##  if (-sum(rate.eigen$stationary*diag(rate.eigen$rate.matrix))*interval.len < 10^-5){
  ##    prob.matrix = diag(rep(1,length(rate.eigen$stationary))) +
  ##      rate.eigen$rate.matrix*interval.len
  ##    }else{
  prob.matrix = rate.eigen$vectors%*%
    diag(exp(rate.eigen$values*interval.len))%*%
    rate.eigen$invvectors
  ##  }
  
  return(prob.matrix)
}

rescale.mc = function(rate.mc, scaling.factor){
  if (!("revmc" %in% class(rate.mc)))
    stop("Error: object \"rate.mc\" is not of class \"revmc\"")
  
  rate.mc$rate.matrix = scaling.factor*rate.mc$rate.matrix
  
  if ("eigen" %in% class(rate.mc))
    rate.mc$values = scaling.factor*rate.mc$values
  
  return(rate.mc)
}

stat.marg.mean.markov.jumps = function(rate.model, regist.matrix, interval.len){
  if (!("revmc" %in% class(rate.model)))
    stop("Error: object \"rate.model\" is not of class \"revmc\"")
  
  one.vector = rep(1, length(rate.model$stationary))
  return(t(rate.model$stationary)%*%(rate.model$rate.matrix*regist.matrix)%*%one.vector*interval.len)
}

joint.mean.markov.jumps = function(rate.eigen, regist.matrix, interval.len){
  
  if (!("eigen" %in% class(rate.eigen)))
    stop("Error: object \"rate.eigen\" is not of class \"eigen\"")
  
  if (!("revmc" %in% class(rate.eigen)))
    stop("Error: object \"rate.eigen\" is not of class \"revmc\"")
  
  if (!is.matrix(regist.matrix) || (nrow(regist.matrix) != ncol(regist.matrix))) 
    stop("\"regist.matrix\" must be a square matrix")
  
  if (nrow(rate.eigen$rate.matrix) != nrow(regist.matrix)) 
    stop("dimensions of rate and register matrices do not match")
  
  if (prod((regist.matrix == 1) + (regist.matrix == 0)) == 0)
    stop("all entries of \"regist.matrix\" must be either 0 or 1") 
  
  space.size = dim(regist.matrix)[1]
  
  zero.matrix = matrix(0,space.size,space.size)
  factorial.moments = zero.matrix
  rate.reg = rate.eigen$rate.matrix*regist.matrix
  
  if (-sum(rate.eigen$stationary*diag(rate.eigen$rate.matrix))*interval.len < 0.001){## if time interval is very small do first order Taylor expansion
    factorial.moments = rate.reg*interval.len +
      (rate.eigen$rate.matrix%*%rate.reg + rate.reg%*%rate.eigen$rate.matrix)*(interval.len)^2/2 +
      rate.eigen$rate.matrix%*%rate.reg%*%rate.eigen$rate.matrix*(interval.len)^3/6
  }else{
    int.matrix = .Call("aux_mat1", rate.eigen$values, interval.len)
    
    factorial.moments = rate.eigen$vectors%*%
      (int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors))%*%
      rate.eigen$invvectors
  }
  
  return(factorial.moments)
}

ind.codon.eigen = function(codon1.eigen, codon2.eigen, codon3.eigen){
  if (!(("revmc" %in% class(codon1.eigen)) && ("eigen" %in% class(codon1.eigen))))
    stop("Error: object \"codon1.mc\" is not of class \"revmc\" and \"eigen\"")
  
  if (!(("revmc" %in% class(codon2.eigen)) && ("eigen" %in% class(codon2.eigen))))
    stop("Error: object \"codon2.mc\" is not of class \"revmc\" and \"eigen\"")
  
  if (!(("revmc" %in% class(codon3.eigen)) && ("eigen" %in% class(codon3.eigen))))
    stop("Error: object \"codon3.mc\" is not of class \"revmc\" and \"eigen\"")
  
  
  rate.mat = kronecker.sum(kronecker.sum(codon1.eigen$rate.matrix, codon2.eigen$rate.matrix),codon3.eigen$rate.matrix)
  
  stat = codon1.eigen$stationary%x%codon2.eigen$stationary%x%codon3.eigen$stationary
  
  ident.vec = rep(1,length(codon1.eigen$stationary))
  
  eigen.val = (codon1.eigen$values%x%ident.vec + ident.vec%x%codon2.eigen$values)%x%
    ident.vec + ident.vec%x%ident.vec%x%codon3.eigen$values
  
  right.eigen.vec = (codon1.eigen$vectors%x%codon2.eigen$vectors)%x%codon3.eigen$vectors
  
  left.eigen.vec = t((t(codon1.eigen$invvectors)%x%t(codon2.eigen$invvectors))%x%
                       t(codon3.eigen$invvectors))
  
  
  markov.chain = list()
  markov.chain[["rate.matrix"]] = rate.mat
  markov.chain[["stationary"]] = stat
  markov.chain[["values"]] = eigen.val
  markov.chain[["vectors"]] = right.eigen.vec
  markov.chain[["invvectors"]] = left.eigen.vec
  class(markov.chain) = c("revmc", "eigen")
  
  return(markov.chain)
  
}


