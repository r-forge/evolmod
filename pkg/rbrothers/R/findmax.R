
findmax<-function(probs,i){
 answer<-2
 for(j in 2:dim(probs)[2]){
  if(probs[i,j]>probs[i,answer]){answer<-j}
 }
 return(answer)
}

