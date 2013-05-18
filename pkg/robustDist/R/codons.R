# 1="A" 2="G" 3="C" 4="T"

TTV<-function(char1,char2)
{
	return(  1*((char1 < 3)&(char2 < 3)) 
		  + 1*((char1 > 2)&(char2 > 2))
		  - 1*((char1 > 2)&(char2 < 3)) 
		  - 1*((char2 > 2)&(char1 < 3)) );
}

# Returns 1 if transistional difference, -1 if transversional, 0 if 
# no substitutions or two or more substitutions

codonTest<-function(codon1, codon2)
{
	if (((codon1[1] != codon2[1])&(codon1[2] != codon2[2]))
		|((codon1[1] != codon2[1])&(codon1[3] != codon2[3]))
		|((codon1[2] != codon2[2])&(codon1[3] != codon2[3]))
		|((codon1[1] == codon2[1])&(codon1[2] == codon2[2])&(codon1[3] == codon2[3])))
	{
		return (0);
	}
	else
	{
		return ((codon1[1] != codon2[1])*TTV(codon1[1],codon2[1])
			+ (codon1[2] != codon2[2])*TTV(codon1[2],codon2[2])
			+ (codon1[3] != codon2[3])*TTV(codon1[3],codon2[3]))
	}	
}

## dna.indices is a vector of size 3 with elements ranging from 0 to 4 

dna.to.codon = function(dna.indices){
  
  if (!is.vector(dna.indices) || (length(dna.indices) != 3))
    stop("\"dna.indices\" must be a vector of length 3")
  
  if((dna.indices[1] == 0)|(dna.indices[2] == 0)|(dna.indices[3] == 0)){
    return(0)
  }else{
    p = (as.numeric(dna.indices[1])-1)*16+(as.numeric(dna.indices[2])-1)*4+(as.numeric(dna.indices[3])-1)
    p = p + 1
    return(p)   
  }	
}


codon.to.dna = function(codon.index){	
  codon.index = codon.index-1;
  
  p = c(0,0,0);
  
  p[1] = floor(codon.index/16)+1;
  codon.index = codon.index - floor(codon.index/16)*16
  p[2] = floor(codon.index/4)+1
  codon.index = codon.index - floor(codon.index/4)*4
  p[3] = floor(codon.index)+1
  
  return(p)
}

is.single.mutation = function(codon.index1, codon.index2){
  
  return.value = FALSE

  dna.index1 = codon.to.dna(codon.index1)
  dna.index2 = codon.to.dna(codon.index2)
  ind.diff = dna.index1 - dna.index2
  
  if (length(ind.diff[ind.diff == 0]) == 2){
    return.value = TRUE
  }

  return(return.value)
}

is.syn.mutation = function(codon.index1, codon.index2){
  return.value = FALSE
  
  if (is.single.mutation(codon.index1, codon.index2) &
      ## Added condition so that stop to stops are not counted
  	  (codon.to.amino(codon.index1) != 0) &  
      (codon.to.amino(codon.index1) == codon.to.amino(codon.index2))){
    return.value = TRUE
  }
  return(return.value)
}


is.nonsyn.mutation = function(codon.index1, codon.index2){
  return.value = FALSE
  
  if (is.single.mutation(codon.index1, codon.index2) &
      ## Added condition so that stop to stops are not counted
  	  (codon.to.amino(codon.index1) != 0) & (codon.to.amino(codon.index2) != 0) &
      (codon.to.amino(codon.index1) != codon.to.amino(codon.index2))){
    return.value = TRUE
  }
  return(return.value)
}

## construct a 64x64 register matrix for synonymous mutations

regist.synonym = function(){
  regist.matrix = matrix(0, nrow = 64, ncol = 64)

  for (i in 1:64){
    for (j in 1:64){
      if (is.syn.mutation(i,j))
        regist.matrix[i,j] = 1
    }
  }
  return(regist.matrix)
}

## construct a 64x64 register matrix for nonsynonymous mutations

regist.nonsynonym = function(){
  regist.matrix = matrix(0, nrow = 64, ncol = 64)

  for (i in 1:64){
    for (j in 1:64){
      if (is.nonsyn.mutation(i,j))
        regist.matrix[i,j] = 1
    }
  }
  return(regist.matrix)
}

## construct a 64x64 register matrix for synonymous transitions

regist.synonym.ts = function(){
  regist.matrix = matrix(0, nrow = 64, ncol = 64)

  for (i in 1:64){
    for (j in 1:64){
      if (is.syn.mutation(i,j) & codonTest(codon.to.dna(i), codon.to.dna(j))==1)
        regist.matrix[i,j] = 1
    }
  }
  return(regist.matrix)
}

## construct a 64x64 register matrix for synonymous transversions

regist.synonym.tv = function(){
  regist.matrix = matrix(0, nrow = 64, ncol = 64)

  for (i in 1:64){
    for (j in 1:64){
      if (is.syn.mutation(i,j) & codonTest(codon.to.dna(i), codon.to.dna(j))==-1)
        regist.matrix[i,j] = 1
    }
  }
  return(regist.matrix)
}

## construct a 64x64 register matrix for nonsynonymous transitions

regist.nonsynonym.ts = function(){
  regist.matrix = matrix(0, nrow = 64, ncol = 64)

  for (i in 1:64){
    for (j in 1:64){
      if (is.nonsyn.mutation(i,j) & codonTest(codon.to.dna(i), codon.to.dna(j))==1)
        regist.matrix[i,j] = 1
    }
  }
  return(regist.matrix)
}

## construct a 64x64 register matrix for nonsynonymous transversions

regist.nonsynonym.tv = function(){
  regist.matrix = matrix(0, nrow = 64, ncol = 64)

  for (i in 1:64){
    for (j in 1:64){
      if (is.nonsyn.mutation(i,j) & codonTest(codon.to.dna(i), codon.to.dna(j))==-1)
        regist.matrix[i,j] = 1
    }
  }
  return(regist.matrix)
}

## Ala = 1
## Arg = 2
## Asn = 3
## Asp = 4
## Cys = 5
## Glu = 6 (E)
## Gln = 7 (Q)
## Gly = 8
## His = 9
## lle = 10
## Leu = 11
## Lys = 12
## Met = 13
## Phe = 14
## Pro = 15
## Ser = 16
## Thr = 17
## Trp = 18
## Tyr = 19
## Val = 20
## STOP = 0

codon.to.amino = function(codon.index){
  amino.acid = switch(codon.index,
    12,
    12,
    3,
    3,
    2,
    2,
    16,
    16,
    17,
    17,
    17,
    17,
    10,
    13, ## start codon
    10,
    10,
    7,    ## Was previously coded as 6 incorrectly
    7,    ## Was previously coded as 6
    4,
    4,
    8,
    8,
    8,
    8,
    1,
    1,
    1,
    1,
    20,
    20,
    20,
    20,
    6,
    6,
    9,
    9,
    2,
    2,
    2,
    2,
    15,
    15,
    15,
    15,
    11,
    11,
    11,
    11,
    0, ## stop codon
    0, ## stop codon
    19,
    19,
    0, ## stop codon
    18,
    5,
    5,
    16,
    16,
    16,
    16,
    11,
    11,
    14,
    14
    )
  
  return(amino.acid)
}


# 14, 49 and 53 correspond to the assignations for the two stop codons 
# and start codons 

codonHKY.model <- function(transition.rate, transversion.rate, stationary, scale = F)
{
	if (length(stationary) !=4)
	{	stop("Problem")	}

	q.matrix = matrix(rep(0,61*61),nrow = 61, ncol = 61);
	
	I = 0;

	for (i in c(1:64))
	{
		if ((i == 14)|(i == 49)|(i == 53))
		{ 
			I = I + 1;
		}
		else
		{
			codonI = convert10to4(i);
			J = 0;

			for (j in c(1:64))
			{
				if ((j == 14)|(j == 49)|(j == 53))
				{
					J = J + 1;
				}
				else
				{
					codonJ = convert10to4(j); 
					test = codonTest(codonI,codonJ);

					if (test != 0)
					{	
						if (test == 1)
						{	q.matrix[i-I,j-J] = transition.rate; }
						else 
						{	q.matrix[i-I,j-J] = transversion.rate;   }
						
						q.matrix[i-I,j-J] = q.matrix[i-I,j-J]*(stationary[codonJ[1]]*(codonI[1] != codonJ[1]) + stationary[codonJ[2]]*(codonI[2] != codonJ[2]) + stationary[codonJ[3]]*(codonI[3] != codonJ[3]))
						}
					else
					{
						q.matrix[i-I,j-J] = 0;
					}
				}
			}
			
		}
	}
	
	 for (i in c(1:61)){
   	 q.matrix[i,i] = -sum(q.matrix[i,])
  }
	return(q.matrix);
	

}

