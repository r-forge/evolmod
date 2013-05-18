as.revmc = function(q.matrix, pi, par.vector=NULL){
  
  if (!is.matrix(q.matrix) || (nrow(q.matrix) != ncol(q.matrix))) 
    stop("\"q.matrix\" must be a square matrix")
  
  if (nrow(q.matrix) != length(pi))
    stop("dimenstions of \"q.matrix\" and \"pi\" do not match")
  
  markov.chain = list()
  markov.chain[["rate.matrix"]] = q.matrix
  markov.chain[["stationary"]] = pi
  markov.chain[["par"]] = par.vector
  
  if (sum((markov.chain$stationary%*%markov.chain$rate.matrix)^2) > 10^-3)
    stop("rate matrix and stationary distribution must satisfy pi%*%q.matrix = 0")
  
  temp.matrix = diag(markov.chain$stationary)%*%markov.chain$rate.matrix
  if (sum((temp.matrix - t(temp.matrix))^2) > 10^-7)
    stop("model is not reversible")
  
  class(markov.chain) = "revmc"
  return(markov.chain)
}


gtr.mc = function(gtr.rates, stationary, scale = F){
  if(length(gtr.rates) != 6)
    stop("\"gtr.rates\" must have length 6")
  
  
  if(length(stationary) != 4)
    stop("vector \"stationary\" must be of length 4")
  
  q.matrix = matrix(c(0,gtr.rates[1], gtr.rates[2], gtr.rates[3], gtr.rates[1], 0, gtr.rates[4],
                      gtr.rates[5], gtr.rates[2],gtr.rates[4],0,gtr.rates[6], gtr.rates[3],gtr.rates[5],gtr.rates[6],0),
                    byrow = T, nrow = 4, ncol = 4)
  
  q.matrix = q.matrix%*%diag(stationary)
  
  for (i in c(1:4)){
    q.matrix[i,i] = -sum(q.matrix[i,])
  }
  
  factor = 1.0
  
  if(scale){
    factor = sum(-diag(q.matrix)*stationary)
  }
  
  return(as.revmc(q.matrix/factor, stationary, gtr.rates))
  
}

f84.mc = function(f84.rates, stationary, scale = F){
  if(length(f84.rates) != 2)
    stop("\"f84.rates\" must have length 2")
  
  ts.rate = f84.rates[1]
  tv.rate = f84.rates[2]
  
  if(length(stationary) != 4)
    stop("vector \"stationary\" must be of length 4")
  
  pi.R = stationary[1] + stationary[2]
  pi.Y = stationary[3] + stationary[4]
  
  q.matrix = matrix(c(0,ts.rate/pi.R+tv.rate, tv.rate, tv.rate, ts.rate/pi.R+tv.rate, 0, tv.rate,
                      tv.rate,tv.rate,tv.rate,0,ts.rate/pi.Y+tv.rate, tv.rate,tv.rate,ts.rate/pi.Y+tv.rate,0),
                    byrow = T, nrow = 4, ncol = 4)
  
  q.matrix = q.matrix%*%diag(stationary)
  
  for (i in c(1:4)){
    q.matrix[i,i] = -sum(q.matrix[i,])
  }
  
  factor = 1.0
  
  if(scale){
    factor = sum(-diag(q.matrix)*stationary)
  }
  
  return(as.revmc(q.matrix/factor, stationary, f84.rates))
}

jc.mc = function(alpha, mc.stat,scale = F){
  my.jc = hky.mc(c(alpha, alpha), mc.stat, scale)
  my.jc$par = alpha
  return(my.jc)
}

hky.mc = function(hky.rates, stationary, scale = F){
  if(length(hky.rates) != 2)
    stop("\"hky.rates\" must have length 2")
  
  ts.rate = hky.rates[1]
  tv.rate = hky.rates[2]
  
  if(length(stationary) != 4)
    stop("vector \"stationary\" must be of length 4")
  
  q.matrix = matrix(c(0,ts.rate, tv.rate, tv.rate, ts.rate, 0, tv.rate,
                      tv.rate,tv.rate,tv.rate,0,ts.rate, tv.rate,tv.rate,ts.rate,0),
                    byrow = T, nrow = 4, ncol = 4)
  
  q.matrix = q.matrix%*%diag(stationary)
  
  for (i in c(1:4)){
    q.matrix[i,i] = -sum(q.matrix[i,])
  }
  
  factor = 1.0
  
  if(scale){
    factor = sum(-diag(q.matrix)*stationary)
  }
  
  return(as.revmc(q.matrix/factor, stationary, hky.rates))
  
}





