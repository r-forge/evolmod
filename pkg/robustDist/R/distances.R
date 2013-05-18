site.pair.likelihood = function(seq.pair, transition.prob){
  return.val = 0
  
  if ((seq.pair[1] != 0) && (seq.pair[2] != 0)){ 
    return.val = transition.prob[seq.pair[1],seq.pair[2]]
  }else{
    return.val = 1.0
  }
  return(return.val)
}

site.pair.counts = function(seq.pair, cond.counts, marg.counts, stat.counts){
  return.val = 0
  
  if ((seq.pair[1] != 0) && (seq.pair[2] != 0)){ 
    return.val = cond.counts[seq.pair[1],seq.pair[2]]
  }else{
    if ((seq.pair[1] == 0) && (seq.pair[2] != 0)){
      return.val = marg.counts[seq.pair[2]]
    }else{
      if ((seq.pair[1] != 0) && (seq.pair[2] == 0)){
        return.val = marg.counts[seq.pair[1]]
      }else{
        return.val = stat.counts
      }
    }
  }
  
  return(return.val)
}


pair.likelihood.hky = function(para,seq.table,pi){

  if(dim(seq.table)[1] != 2)
    stop("number of rows in \"seq.table\" must be 2")

  seqLen = dim(seq.table)[2]
  
  rev.model = hky.model(para[1],para[2], pi, scale = F)	
  
  rev.eigen = as.eigen(rev.model)
  
  rev.prob = matexp(rev.eigen,1)
  
  like = apply(seq.table, 2, FUN = site.pair.likelihood,
    transition.prob = rev.prob)

  
  return (sum(log(like)))
}

pair.log.like.old = function(mc.rates, mc.stat, mc.model, seq.table){
  
  if(dim(seq.table)[1] != 2)
    stop("number of rows in \"seq.table\" must be 2")
  
  seqLen = dim(seq.table)[2]
  
  rev.model = mc.model(mc.rates, mc.stat, scale=F)	
  
  rev.eigen = as.eigen(rev.model)
  
  rev.prob = matexp(rev.eigen,1)
  
  like = apply(seq.table, 2, FUN = site.pair.likelihood,
    transition.prob = rev.prob)
  
  return (sum(log(like)))
}

pair.log.like = function(mc.log.rates, mc.stat, mc.model, site.patterns){
  
  if(dim(site.patterns)[1] != length(mc.stat))
    stop("number of rows and columns in \"site.patterns\" must be length of mc.stat")
  
  rev.model = mc.model(exp(mc.log.rates), mc.stat, scale=F)	
  
  rev.eigen = as.eigen(rev.model)
  
  rev.prob = matexp(rev.eigen,1)
  
  like = sum(log(rev.prob)*site.patterns)
  
  return(like)
}



## 1 = A, 2 = G, 3 = C, 4 = T

fit.hky = function(mc.stat, seq.table){

  hky.eigen = as.eigen.hky(mc.stat,1,1)

  pair.counts = count.pair.patterns(seq.table,4)

  pi.Y = mc.stat[3] + mc.stat[4]
  pi.R = mc.stat[1] + mc.stat[2]

  return(c(aux.sum2,aux.sum3,aux.sum4))
}  


fit.f84 = function(mc.stat, seq.table){

  pi.R = mc.stat[1] + mc.stat[2]
  pi.Y = mc.stat[3] + mc.stat[4]

  ## compute the frequenies of transitions and tranversions

  pair.counts = .Call("pair_patterns",seq.table,4)
  ##pair.counts = count.pair.patterns(seq.table,4)

  seq.length = dim(seq.table)[2]
  
  ts.freq = (pair.counts[1,2] + pair.counts[2,1] + pair.counts[3,4] + pair.counts[4,3])/seq.length
  tv.freq = (pair.counts[1,3] + pair.counts[1,4] + pair.counts[2,3] + pair.counts[2,4] +
             pair.counts[3,1] + pair.counts[3,2] + pair.counts[4,1] + pair.counts[4,2])/seq.length
  
  return.object = NULL

  if (ts.freq*tv.freq == 0){
    return.object = jc.mc(-log(1-0.75*(ts.freq+tv.freq)), rep(0.25,4),scale = F)
  }else{
  
    A = mc.stat[1]*mc.stat[2]/pi.R + mc.stat[3]*mc.stat[4]/pi.Y
    B = mc.stat[1]*mc.stat[2] + mc.stat[3]*mc.stat[4]
    C = pi.R*pi.Y
    
    
    tv.rate = -log(1-tv.freq/(2*C))
    ts.rate = -log(1 - B/A - exp(tv.rate)*(0.5*ts.freq-B)/A)
    
    return.object = f84.mc(c(ts.rate,tv.rate),mc.stat,F)
  }
  
  return(return.object)
}

fit.poisson = function(seq.table){

  ## compute the frequenies of transitions and tranversions

  pair.counts = count.pair.patterns(seq.table,2)

  seq.length = dim(seq.table)[2]
  
  invariable.freq = (pair.counts[1,1] + pair.counts[2,2])/seq.length
  variable.freq = 1 - invariable.freq

  return.object = NULL

  my.rate = -log(invariable.freq-variable.freq)
  
  return.object = two.state.mc(my.rate,c(0.5,0.5),F)

  return(return.object)
}

fit.jc = function(seq.table){

  ## compute the frequenies of transitions and tranversions

  pair.counts = .Call("pair_patterns",seq.table,4)##count.pair.patterns(seq.table,4)

  seq.length = dim(seq.table)[2]
  
  variable.freq = 1-sum(diag(pair.counts))/seq.length

  return.object = NULL

  my.rate = -log(1-4/3*variable.freq)
  
  return.object = jc.mc(my.rate,c(0.25,0.25,0.25,0.25),F)
  
  return(return.object)
}

get.f84.rates = function(mc.stat, seq.table){

  pi.R = mc.stat[1] + mc.stat[2]
  pi.Y = mc.stat[3] + mc.stat[4]

  ## compute the frequenies of transitions and tranversions

  pair.counts = count.pair.patterns(seq.table,4)

  seq.length = dim(seq.table)[2]
  
  ts.freq = (pair.counts[1,2] + pair.counts[2,1] + pair.counts[3,4] + pair.counts[4,3])/seq.length
  tv.freq = (pair.counts[1,3] + pair.counts[1,4] + pair.counts[2,3] + pair.counts[2,4] +
             pair.counts[3,1] + pair.counts[3,2] + pair.counts[4,1] + pair.counts[4,2])/seq.length

  A = mc.stat[1]*mc.stat[2]/pi.R + mc.stat[3]*mc.stat[4]/pi.Y
  B = mc.stat[1]*mc.stat[2] + mc.stat[3]*mc.stat[4]
  C = pi.R*pi.Y
    
    
  tv.rate = -log(1-tv.freq/(2*C))
  ts.rate = -log(1 - B/A - exp(tv.rate)*(0.5*ts.freq-B)/A)

  return(c(ts.rate,tv.rate))
  
}


fit.mc = function(mc.model, mc.rate.init, mc.stat, seq.table){
  return.val = NULL

  ## determine counts of nucleotide patters,
  ## can ignore patters with gaps as they do not contribute to the likelihood
  
  optim.out = optim(mc.rate.init, pair.log.like,
    lower = rep(-10,length(mc.rate.init)), upper = rep(10,length(mc.rate.init)),
    mc.stat = mc.stat, mc.model = mc.model,
    site.patterns = count.pair.patterns(seq.table,length(mc.stat)),
    control = list(fnscale = -1), method = "L-BFGS-B")
  
  ##print(optim.out)
  
  if (optim.out$convergence != 0){
    return.val = FALSE
  }else{
    return.val = mc.model(exp(
      optim.out$par
    ), mc.stat, scale = F)
  }

  return(return.val)
}


count.pair.patterns = function(seq.table, space.size){
  
  pattern.matrix = diag(rep(0,space.size))
  
  for (i in 1:space.size){
    for (j in 1:space.size){
      pattern.matrix[i,j] = sum((seq.table[1,]==i)*(seq.table[2,]==j))
    }
  }

  return(pattern.matrix)
}

count.pair.patterns.old = function(seq.table, space.size){
  
  pattern.matrix = diag(rep(0,space.size))
  
  for (i in 1:dim(seq.table)[2]){
    if ((seq.table[1,i] != 0) && (seq.table[2,i] != 0)){
      pattern.matrix[seq.table[1,i],seq.table[2,i]] = pattern.matrix[seq.table[1,i],seq.table[2,i]] + 1
    }
  }

  return(pattern.matrix)
}

extended.count.pair.patterns = function(seq.table, space.size){
  
  pattern.matrix = diag(rep(0,space.size+1))
  
  for (i in 0:space.size){
    for (j in 0:space.size){
      pattern.matrix[i+1,j+1] = sum((seq.table[1,]==i)*(seq.table[2,]==j))
    }
  }

  return(pattern.matrix)
}

pair.conv.dist = function(rev.model, regist.matrix){
  return(as.numeric(stat.marg.mean.markov.jumps(rev.model, regist.matrix, 1.0)))
}

pair.robust.dist = function(rev.eigen, regist.matrix, seq.table){

  space.size = dim(regist.matrix)[1]
  seq.length = dim(seq.table)[2]

  joint.counts = joint.mean.markov.jumps(rev.eigen, regist.matrix, 1.0)
  rev.prob = matexp(rev.eigen,1)
  cond.counts = joint.counts/rev.prob
  
  one.vector = rep(1,nrow(joint.counts))
  marg.counts = as.vector(joint.counts%*%one.vector)
  stat.counts = as.numeric(t(rev.eigen$stationary)%*%marg.counts)

  ##site.pat = .Call("ext_pair_patterns",seq.table,space.size)
  ##site.pat = extended.count.pair.patterns(seq.table,space.size)

  ##total.counts = sum(site.pat[2:(space.size+1),2:(space.size+1)]*cond.counts) +
  ##  sum(site.pat[1,2:(space.size+1)]*marg.counts) + sum(site.pat[2:(space.size+1),1]*marg.counts) +
  ##    stat.counts*site.pat[1,1]

  total.counts = .Call("pair_counts",seq.table,cond.counts,marg.counts,stat.counts)
  
  return(total.counts/seq.length)
}

pair.robust.dist.old = function(rev.eigen, regist.matrix, seq.table){

  joint.counts = joint.mean.markov.jumps(rev.eigen, regist.matrix, 1.0)
  rev.prob = matexp(rev.eigen,1)
  cond.counts = joint.counts/rev.prob
  
  one.vector = rep(1,nrow(joint.counts))
  marg.counts = joint.counts%*%one.vector
  stat.counts = as.numeric(t(rev.eigen$stationary)%*%marg.counts)
  
  
  site.counts = apply(seq.table, 2, FUN = site.pair.counts,
    cond.counts = cond.counts, marg.counts = marg.counts, stat.counts = stat.counts)
  
  return(mean(site.counts))
}

conv.dist = function(mc.model, mc.rate.init, mc.size, regist.matrix, seq.table){

  ## estimate stationary distribution
  mc.stat = rep(0,mc.size)
  
  for (i in 1:mc.size){
    ## count occurences of i
    mc.stat[i] = length(which(seq.table == i))
  }

  mc.stat = mc.stat/sum(mc.stat)

  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){
      ## bind sequences i and j into a matrix
      pair.seq = rbind(seq.table[i,],seq.table[j,])

      ## fit markov chain model using sequences i and j

      pair.mc = fit.mc(mc.model, mc.rate.init, mc.stat, pair.seq)

      ## compute distance between i and j using the fitted model

      dist.matrix[j,i] = pair.conv.dist(pair.mc, regist.matrix)
    }
  }

  return(dist.matrix)
}

conv.f84 = function(regist.matrix, seq.table){

  mc.size = 4
  
  ## estimate stationary distribution
  mc.stat = emp.freq(seq.table,mc.size)

  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){
      ## bind sequences i and j into a matrix
      pair.seq = rbind(seq.table[i,],seq.table[j,])

      ## fit markov chain model using sequences i and j

      pair.mc = fit.f84(mc.stat, pair.seq)

      ## compute distance between i and j using the fitted model
      dist.matrix[j,i] = pair.conv.dist(pair.mc, regist.matrix)
      dist.matrix[i,j] = dist.matrix[j,i]
    }
  }

  return(dist.matrix)
}

conv.jc = function(regist.matrix, seq.table){

  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){
      ## bind sequences i and j into a matrix
      pair.seq = rbind(seq.table[i,],seq.table[j,])

      ## fit markov chain model using sequences i and j

      pair.mc = fit.jc(pair.seq)

      ## compute distance between i and j using the fitted model
      dist.matrix[j,i] = pair.conv.dist(pair.mc, regist.matrix)
      dist.matrix[i,j] = dist.matrix[j,i]
    }
  }

  return(dist.matrix)
}


robust.f84 = function(regist.matrix, seq.table){

  ## estimate stationary distribution
  mc.stat = emp.freq(seq.table,mc.size)

  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){
      ## bind sequences i and j into a matrix
      pair.seq = rbind(seq.table[i,],seq.table[j,])

      ## fit markov chain model using sequences i and j

      pair.mc = fit.f84(mc.stat, pair.seq)

      ## compute eigen decomposition of the Markov generator

      pair.mc.eigen = as.eigen(pair.mc)
      
      ## compute distance between i and j using the fitted model
      
      dist.matrix[j,i] = pair.robust.dist(pair.mc.eigen, regist.matrix, pair.seq)
      dist.matrix[i,j] = dist.matrix[j,i]
    }
  }
  
  return(dist.matrix)
}

robust.jc = function(regist.matrix, seq.table){
  
  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){
      ## bind sequences i and j into a matrix
      pair.seq = rbind(seq.table[i,],seq.table[j,])

      ## fit markov chain model using sequences i and j

      pair.mc = fit.jc(mc.stat, pair.seq)

      ## compute eigen decomposition of the Markov generator
      pair.mc.eigen = as.eigen(pair.mc)
      
      ## compute distance between i and j using the fitted model
      
      dist.matrix[j,i] = pair.robust.dist(pair.mc.eigen, regist.matrix, pair.seq)
      dist.matrix[i,j] = dist.matrix[j,i]
    }
  }
  
  return(dist.matrix)
}


robust.dist = function(mc.model, mc.rate.init, mc.size, regist.matrix, seq.table){
  
  ## estimate stationary distribution
  mc.stat = rep(0,mc.size)
  
  for (i in 1:mc.size){
    ## count occurences of i
    mc.stat[i] = length(which(seq.table == i))
  }

  mc.stat = mc.stat/sum(mc.stat)

  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){
      ## bind sequences i and j into a matrix
      pair.seq = rbind(seq.table[i,],seq.table[j,])

      ## fit markov chain model using sequences i and j

      pair.mc = fit.mc(mc.model, mc.rate.init, mc.stat, pair.seq)

      ## compute eigen decomposition of the Markov generator

      pair.mc.eigen = as.eigen(pair.mc)
      
      ## compute distance between i and j using the fitted model
      
      dist.matrix[j,i] = pair.robust.dist(pair.mc.eigen, regist.matrix, pair.seq)
    }
  }
  
  return(dist.matrix)
}


## regist.matrix is on the codon space

robust.codon.dist = function(mc.model, mc.rate.init, mc.size, regist.matrix, seq.table){

  dna.align.len = dim(seq.table)[2]
  
  ## make a codon table 

  codon.align.obj = dna.to.codon.align(seq.table)
  
  ## estimate stationary distribution
  mc.stat = rep(0,mc.size)
  
  for (i in 1:mc.size){
    ## count occurences of i
    mc.stat[i] = length(which(seq.table == i))
  }

  mc.stat = mc.stat/sum(mc.stat)

  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)

  ## make individual codon alignments

  seq.table.codon1 = seq.table[,seq(from = 1, to = dna.align.len-2, by = 3)]
  seq.table.codon2 = seq.table[,seq(from = 2, to = dna.align.len-1, by = 3)]
  seq.table.codon3 = seq.table[,seq(from = 3, to = dna.align.len, by = 3)]
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){

      
      ## bind sequences i and j into a matrix of codon position 1
      pair.seq.codon1 = rbind(seq.table.codon1[i,],seq.table.codon1[j,])
      ## bind sequences i and j into a matrix of codon position 2
      pair.seq.codon2 = rbind(seq.table.codon2[i,],seq.table.codon2[j,])
      ## bind sequences i and j into a matrix of codon position 3
      pair.seq.codon3 = rbind(seq.table.codon3[i,],seq.table.codon3[j,])

      pair.seq = rbind(codon.align.obj[i,],codon.align.obj[j,])
      
      ## fit 3 markov chain models using sequences i and j

      pair.mc.codon1 = fit.mc(mc.model, mc.rate.init, mc.stat, pair.seq.codon1)
      pair.mc.codon2 = fit.mc(mc.model, mc.rate.init, mc.stat, pair.seq.codon2)
      pair.mc.codon3 = fit.mc(mc.model, mc.rate.init, mc.stat, pair.seq.codon3)


      ## compute eigen decomposition of the three Markov generators

      ##pair.mc.eigen = as.eigen(pair.mc)

      ## make codon matrix from 3 individual dna generators

      pair.codon = ind.codon.mc( pair.mc.codon1,  pair.mc.codon2,  pair.mc.codon3)

      pair.codon.eigen = as.eigen(pair.codon)
      
      ## compute distance between i and j using the fitted model
      
      dist.matrix[j,i] = pair.robust.dist(pair.codon.eigen, regist.matrix, pair.seq)
      dist.matrix[i,j] = dist.matrix[j,i]
    }
  }
  
  return(dist.matrix)
}


## regist.matrix is on the codon space

robust.codon.dist.new = function(mc.model, mc.rate.init, mc.size, regist.matrix, seq.table){

  dna.align.len = dim(seq.table)[2]

   ## make a codon table 

  codon.align.obj = dna.to.codon.align(seq.table)
  
  ## estimate stationary distribution
  mc.stat = rep(0,mc.size)
  
  for (i in 1:mc.size){
    ## count occurences of i
    mc.stat[i] = length(which(seq.table == i))
  }

  mc.stat = mc.stat/sum(mc.stat)

  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)

  convergence = TRUE

  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){

      pair.dna = rbind(seq.table[i,],seq.table[j,])
      pair.codon = rbind(codon.align.obj[i,],codon.align.obj[j,])
      
      ## fit 3 markov chain models using sequences i and j

      pair.mc.codon = fit.mc(mc.model, mc.rate.init, mc.stat, pair.dna)

      
      if (!("revmc" %in% class(pair.mc.codon))){
        convergence = FALSE
      }
      
      ## make codon matrix from 3 individual dna generators

      pair.mc.eigen = as.eigen(pair.mc.codon)
      
      pair.codon.eigen = ind.codon.eigen(pair.mc.eigen,  pair.mc.eigen,  pair.mc.eigen)

      ## compute distance between i and j using the fitted model
      
      dist.matrix[j,i] = pair.robust.dist(pair.codon.eigen, regist.matrix, pair.codon)
      dist.matrix[i,j] = dist.matrix[j,i]
    }
  }

  if (convergence == TRUE){
    return(dist.matrix)
  }else{
    return(FALSE)
  }
}


robust.codon.dist.onef84 = function(seq.table,codon.table,regist.matrix){

  dna.align.len = dim(seq.table)[2]
  
  ## make a codon table 

  codon.align.obj = codon.table
  
  ## estimate stationary distribution
  mc.stat = rep(0,4)
  
  for (i in 1:4){
    ## count occurences of i
    mc.stat[i] = length(which(seq.table == i))
  }

  mc.stat = mc.stat/sum(mc.stat)

  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(0, nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)

  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){

      pair.dna = rbind(seq.table[i,],seq.table[j,])
      pair.codon = rbind(codon.align.obj[i,],codon.align.obj[j,])
      
      ## fit 3 markov chain models using sequences i and j

      pair.mc.codon = fit.f84(mc.stat, pair.dna)
  

      if (prod(pair.mc.codon$rate.matrix==0) == 0){
        
        ## make codon matrix from 3 individual dna generators        
        pair.mc.eigen = as.eigen(pair.mc.codon)
        
        
        pair.codon.eigen = ind.codon.eigen(pair.mc.eigen,  pair.mc.eigen,  pair.mc.eigen)
        
        ## compute distance between i and j using the fitted model        
        dist.matrix[j,i] = pair.robust.dist(pair.codon.eigen, regist.matrix, pair.codon)
        dist.matrix[i,j] = dist.matrix[j,i]
      }
    }
  }

  return(dist.matrix)
}


robust.codon.dist.f84.pos.1.2.3 = function(seq.table, codon.align.obj, regist.matrix){

  dna.align.len = dim(seq.table)[2]
  
  ## make individual codon-position alignments
  
  seq.table.codon1 = seq.table[,seq(from = 1, to = dna.align.len-2, by = 3)]
  seq.table.codon2 = seq.table[,seq(from = 2, to = dna.align.len-1, by = 3)]
  seq.table.codon3 = seq.table[,seq(from = 3, to = dna.align.len, by = 3)]
  
  ## estimate stationary distribution for each codon
  
  mc.stat1 = emp.freq(seq.table.codon1,4)
  mc.stat2 = emp.freq(seq.table.codon2,4)
  mc.stat3 = emp.freq(seq.table.codon3,4)
  
  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)

  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){

      ## bind sequences i and j into a matrix of codon position 1
      pair.seq.codon1 = rbind(seq.table.codon1[i,],seq.table.codon1[j,])
      
      ## bind sequences i and j into a matrix of codon position 2
      pair.seq.codon2 = rbind(seq.table.codon2[i,],seq.table.codon2[j,])
      
      ## bind sequences i and j into a matrix of codon position 3
      pair.seq.codon3 = rbind(seq.table.codon3[i,],seq.table.codon3[j,])

      pair.codon = rbind(codon.align.obj[i,],codon.align.obj[j,])
      
      ## fit 3 markov chain models using sequences i and j

      pair.mc.codon1 = fit.f84(mc.stat1, pair.seq.codon1)
      pair.mc.codon2 = fit.f84(mc.stat2, pair.seq.codon2)
      pair.mc.codon3 = fit.f84(mc.stat3, pair.seq.codon3)
      
      ## make codon matrix from 3 individual dna generators

      pair.mc.eigen1 = as.eigen(pair.mc.codon1)
      pair.mc.eigen2 = as.eigen(pair.mc.codon2)
      pair.mc.eigen3 = as.eigen(pair.mc.codon3)
      
      pair.codon.eigen = ind.codon.eigen(pair.mc.eigen1,  pair.mc.eigen2,  pair.mc.eigen3)

      ## compute distance between i and j using the fitted model
      
      dist.matrix[j,i] = pair.robust.dist(pair.codon.eigen, regist.matrix, pair.codon)
      dist.matrix[i,j] = dist.matrix[j,i]

    }
  }
  
  return(dist.matrix)
}

robust.codon.dist.f84.pos.12.3 = function(seq.table, codon.align.obj, regist.matrix){

  dna.align.len = dim(seq.table)[2]
  
  ## make individual codon-position alignments
  
  seq.table.codon12 = seq.table[,-seq(from = 3, to = dna.align.len, by = 3)]
  seq.table.codon3 = seq.table[,seq(from = 3, to = dna.align.len, by = 3)]
  
  ## estimate stationary distribution for each codon
  
  mc.stat12 = emp.freq(seq.table.codon12,4)
  mc.stat3 = emp.freq(seq.table.codon3,4)
  
  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)

  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){

      ## bind sequences i and j into a matrix of codon positions 1 and 2
      pair.seq.codon12 = rbind(seq.table.codon12[i,],seq.table.codon12[j,])
      
      ## bind sequences i and j into a matrix of codon position 3
      pair.seq.codon3 = rbind(seq.table.codon3[i,],seq.table.codon3[j,])

      pair.codon = rbind(codon.align.obj[i,],codon.align.obj[j,])
      
      ## fit 3 markov chain models using sequences i and j

      pair.mc.codon12 = fit.f84(mc.stat12, pair.seq.codon12)
      pair.mc.codon3 = fit.f84(mc.stat3, pair.seq.codon3)
      
      ## make codon matrix from 3 individual dna generators

      pair.mc.eigen12 = as.eigen(pair.mc.codon12)
      pair.mc.eigen3 = as.eigen(pair.mc.codon3)
      
      pair.codon.eigen = ind.codon.eigen(pair.mc.eigen12,  pair.mc.eigen12,  pair.mc.eigen3)

      ## compute distance between i and j using the fitted model
      
      dist.matrix[j,i] = pair.robust.dist(pair.codon.eigen, regist.matrix, pair.codon)
      dist.matrix[i,j] = dist.matrix[j,i]

    }
  }
  
  return(dist.matrix)
}


robust.codon.dist.f84.pos.123 = function(seq.table, codon.align.obj, regist.matrix){

  dna.align.len = dim(seq.table)[2]
  
  ## estimate stationary distribution
  
  mc.stat = emp.freq(seq.table,4)

  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dist.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  rownames(dist.matrix) = rownames(seq.table)
  colnames(dist.matrix) = rownames(seq.table)
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){

      ## bind sequences i and j into a matrix of codon positions 1 and 2
      pair.seq.dna = rbind(seq.table[i,],seq.table[j,])
      
      pair.codon = rbind(codon.align.obj[i,],codon.align.obj[j,])
      
      ## fit markov chain model using sequences i and j

      pair.mc.dna = fit.f84(mc.stat, pair.seq.dna)
      
      ## make codon matrix from 3 individual dna generators

      pair.mc.eigen = as.eigen(pair.mc.dna)
      
      pair.codon.eigen = ind.codon.eigen(pair.mc.eigen,  pair.mc.eigen,  pair.mc.eigen)

      ## compute distance between i and j using the fitted model
      
      dist.matrix[j,i] = pair.robust.dist(pair.codon.eigen, regist.matrix, pair.codon)
      dist.matrix[i,j] = dist.matrix[j,i]

    }
  }
  
  return(dist.matrix)
}


robust.codon.dist.f84 = function(seq.table, regist.matrix, partition=1, codon.align.obj=NULL){
  return.table = NULL

  ## make a codon table if it is not provided

  if (is.null(codon.align.obj)){
    codon.align.obj = dna.to.codon.align(seq.table)
  }
  
  if (partition == 1){
    return.table = robust.codon.dist.f84.pos.123(seq.table,codon.align.obj,regist.matrix)
  }else{
    if (partition == 2){
      return.table = robust.codon.dist.f84.pos.12.3(seq.table,codon.align.obj, regist.matrix)
    }else{
      return.table = robust.codon.dist.f84.pos.1.2.3(seq.table,codon.align.obj, regist.matrix)
    }
  }

  return(return.table)
}

## regist.matrix is on the codon space

robust.dN.dS = function(seq.table, regist.syn, regist.nonsyn){

  dna.align.len = dim(seq.table)[2]
  
  ## make a codon table 

  codon.align.obj = dna.to.codon.align(seq.table)
  
  ## estimate stationary distribution for each codon
  mc.stat = numeric(4)
  mc.stat1 = numeric(4)
  mc.stat2 = numeric(4)
  mc.stat3 = numeric(4)
  
  for (i in 1:4){
    ## count occurences of i
    mc.stat[i] = length(which(seq.table == i))
    mc.stat1[i] = length(which(seq.table[,seq(from = 1, to = dna.align.len-2, by = 3)] == i))
    mc.stat2[i] = length(which(seq.table[,seq(from = 2, to = dna.align.len-1, by = 3)] == i))
    mc.stat3[i] = length(which(seq.table[,seq(from = 3, to = dna.align.len, by = 3)] == i))
  }

  mc.stat = mc.stat/sum(mc.stat)
  mc.stat1 = mc.stat1/sum(mc.stat1)
  mc.stat2 = mc.stat2/sum(mc.stat2)
  mc.stat3 = mc.stat3/sum(mc.stat3)

  codon.stat = numeric(64)

  for (i in 1:64){
    ## count occurences of i
    codon.stat[i] = length(which(codon.align.obj == i))
  }

  codon.stat = codon.stat/sum(codon.stat)
  
  ## estimate pairwise distances
  seq.num = nrow(seq.table)
  
  dS.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  dN.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  Omega.matrix = matrix(rep(0,seq.num^2), nrow = seq.num, ncol = seq.num)
  
  rownames(dS.matrix) = rownames(seq.table)
  colnames(dS.matrix) = rownames(seq.table)

  rownames(dN.matrix) = rownames(seq.table)
  colnames(dN.matrix) = rownames(seq.table)
  
  rownames(Omega.matrix) = rownames(seq.table)
  colnames(Omega.matrix) = rownames(seq.table)
  
  for (i in 1:(seq.num-1)){
    for (j in (i+1):seq.num){

      pair.dna = rbind(seq.table[i,],seq.table[j,])
      pair.codon = rbind(codon.align.obj[i,],codon.align.obj[j,])
      
      ## fit 1 markov chain model using sequences i and j

      f84.rates = get.f84.rates(mc.stat3, pair.dna[,seq(from = 3, to = dna.align.len, by = 3)])
      
      ##      hky.kappa = codon.hky.like.mc(f84.rates, mc.stat3, codon.stat, scale = F)

      
      ## form 3 markov chains using the same nucleotide rates but
      ## codon-specific stationary distributions
      
      pair.mc.codon1 = f84.mc(f84.rates, mc.stat, F)
      pair.mc.codon2 = f84.mc(f84.rates, mc.stat2, F)
      pair.mc.codon3 = f84.mc(f84.rates, mc.stat3, F)
      
      ## make codon matrix from 3 individual dna generators

      pair.mc.eigen1 = as.eigen(pair.mc.codon1)
      pair.mc.eigen2 = as.eigen(pair.mc.codon2)
      pair.mc.eigen3 = as.eigen(pair.mc.codon3)
      
      
      pair.codon.eigen1 = ind.codon.eigen(pair.mc.eigen1,  pair.mc.eigen1,  pair.mc.eigen1)
      pair.codon.eigen2 = ind.codon.eigen(pair.mc.eigen3,  pair.mc.eigen3,  pair.mc.eigen3)

      ## compute distance between i and j using the fitted model

      S = pair.robust.dist(pair.codon.eigen1, regist.syn, pair.codon)
      N = pair.robust.dist(pair.codon.eigen1, regist.nonsyn, pair.codon)

      seq.dist = (S + N)/3
      
      dS.matrix[j,i] = seq.dist*S/pair.conv.dist(pair.codon.eigen2, regist.syn)
      dS.matrix[i,j] = dS.matrix[j,i]
      
      dN.matrix[j,i] = seq.dist*N/pair.conv.dist(pair.codon.eigen2, regist.nonsyn)
      dN.matrix[i,j] = dN.matrix[j,i]
      
      Omega.matrix[j,i] = dN.matrix[j,i]/dS.matrix[j,i]
      Omega.matrix[i,j] = Omega.matrix[j,i]
      
    }
  }

  return.object = list(dN.matrix, dS.matrix, Omega.matrix)
  names(return.object) = c("dN", "dS", "Omega")
  
  return(return.object)
}

