# This is a wrapper of the "ape" read.dna function. Instead of list of characters for each taxon
# it returns a numerical table with 0 = "-", 1 = "A", 2 = "G", 3 = "C", 4 = "T"

read.phylip = function(filename, format, skip = 0, nlines = 0, comment.char = "#", as.character) {

  if (as.character){
    new.align = read.dna(filename, format, skip = 0, nlines = 0, comment.char = "#", as.character = T)
  }else{
    new.align = read.dna(filename, format, skip = 0, nlines = 0, comment.char = "#")
  }

  if (format == "fasta"){
    num.taxa = length(new.align)
        
    not.in.alphabet = (new.align[[1]] != "a") & (new.align[[1]] != "g")  & (new.align[[1]] != "c") & (new.align[[1]] != "t")

    new.align[[1]][not.in.alphabet] = 0
    new.align[[1]][new.align[[1]] == "a"] = 1
    new.align[[1]][new.align[[1]] == "g"] = 2
    new.align[[1]][new.align[[1]] == "c"] = 3
    new.align[[1]][new.align[[1]] == "t"] = 4
  
    align.table = rbind(as.numeric(new.align[[1]]))
    
    for (i in c(2:num.taxa)){
      not.in.alphabet = (new.align[[i]] != "a") & (new.align[[i]] != "g")  & (new.align[[i]] != "c") & (new.align[[i]] != "t")
      
      new.align[[i]][not.in.alphabet] = 0
      new.align[[i]][new.align[[i]] == "a"] = 1
      new.align[[i]][new.align[[i]] == "g"] = 2
      new.align[[i]][new.align[[i]] == "c"] = 3
      new.align[[i]][new.align[[i]] == "t"] = 4
      
      align.table = rbind(align.table, as.numeric(new.align[[i]]))
    }
    rownames(align.table) = names(new.align)
  }else{

    not.in.alphabet = (new.align != "a")*(new.align != "g")*(new.align != "c")*(new.align != "t")
    not.in.alphabet = not.in.alphabet == 1

    num.row = nrow(new.align)
    num.col = ncol(new.align)

    align.table = matrix(rep(0,num.row*num.col), nrow = num.row, ncol = num.col)
    rownames(align.table) = rownames(new.align)
    ##print(not.in.alphabet)
    
    align.table[not.in.alphabet] = 0
    align.table[new.align == "a"] = 1
    align.table[new.align == "g"] = 2
    align.table[new.align == "c"] = 3
    align.table[new.align == "t"] = 4

  }

  # FIX ME: check if the data were read correctly

  
  return(align.table)  
}


increment<-function(dna){
  
  num.taxa = dim(dna)[1]
  
  inc = rep(1,num.taxa);
  
  for (i in 1:num.taxa){
    
    L = dna[i,]
    
    j = 1;

    while(j < length(L)-2){
      if ((L[j]==1)&(L[j+1]==4)&(L[j+2]==2)){
        inc[i] = j;
        j = length(L);	
      }
      else
        {
          j = j +1;
        }
    }

    if (j == (length(L)-1)){
      inc[i] = "NA"
    }
  }
  
  inc;
}



dna.to.codon.align = function(align.table) {
  
  num.taxa = dim(align.table)[1]
  align.len = dim(align.table)[2]

  if (align.len %% 3 != 0)
    stop("number of columns in \"align.table\" must be multiple of 3") 

  codon.len = align.len/3
  
  codon.table = matrix(0, nrow = num.taxa, ncol = codon.len)
  rownames(codon.table) = rownames(align.table)
  
  for (i in 1:codon.len){
    dna.columns = cbind(align.table[,i*3-2], align.table[,i*3-1], align.table[,i*3])
    codon.table[,i] = apply(dna.columns, MARGIN = 1, dna.to.codon)
  }

  return(codon.table)
}

codon.to.dna.align = function(codon.table) {
  
  num.taxa = dim(codon.table)[1]
  codon.len = dim(codon.table)[2]

  
  dna.len = 3*codon.len
  
  dna.table = matrix(0, nrow = num.taxa, ncol = dna.len)
  rownames(dna.table) = rownames(codon.table)
  
  for (i in 1:codon.len){
    dna.table[,(3*i-2):(3*i)] = t(apply(matrix(codon.table[,i],nrow=num.taxa,ncol=1), MARGIN = 1, codon.to.dna))
  }
  

  return(dna.table)
}



## count occurences of unique sites in the alignment

count.unique.sites = function(my.align){
  
  unique.sites = unique(my.align, MARGIN = 2)
  
  unique.counts = rep(0, dim(unique.sites)[2])
  
  for (i in (1:dim(unique.sites)[2])){
    for (j in (1:dim(my.align)[2])){
      if (sum((my.align[,j] - unique.sites[,i])^2) == 0){
        unique.counts[i] = unique.counts[i] + 1
      }
    }
  }
  
  return(list(unique.sites, unique.counts))
}

## compute empirical frequencies of sequence states

emp.freq = function(seq.table, space.size){
  emp.counts = numeric(space.size)

  for (i in 1:space.size){
    emp.counts[i] = length(which(seq.table == i))
  }

  return(emp.counts/sum(emp.counts))
}

## resample codons with replacement in the nucleotide and codon alignments

resample.codon = function(dna.table, codon.table){

  return.object = list(2)

  dna.length = dim(dna.table)[2]
  
  codon.sample = sample(seq(from=1,to=dna.length, by=3),dna.length/3,replace = T)
  resample.align = dna.table
  resample.codon.align = codon.table
  
  codon1 = seq(from=1,to=dna.length,by=3)
  codon2 = seq(from=2,to=dna.length,by=3)
  codon3 = seq(from=3,to=dna.length,by=3)
  
  resample.align[,codon1] = dna.table[,codon.sample]
  resample.align[,codon2] = dna.table[,codon.sample+1]
  resample.align[,codon3] = dna.table[,codon.sample+2]

  resample.codon.align = codon.table[,(codon.sample+2)/3]

  return.object[[1]] = resample.align
  return.object[[2]] = resample.codon.align

  names(return.object) = c("NucleotideAlignment","CodonAlignment")

  return(return.object)
}
