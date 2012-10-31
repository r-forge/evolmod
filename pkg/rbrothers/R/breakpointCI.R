breakpointCI<-function(x,leftedge,rightedge){

 makebplist<-function(ll,i){
  ii = unlist(strsplit(ll[[i]]," "))
  ii = ii[which(sapply(ii,nchar,USE.NAMES=FALSE)!=0)]
  # tree count:
  tc<-as.numeric(ii[2])
  # potential break points:
  pbp<-as.integer(ii[(5+8*tc):(3+tc+8*tc)])
  # list of trees as strings:
  tl<-ii[4+(0:(tc-1))*8]
  # actual break points
  breaks<-NULL
  for(i in 2:length(tl)) if(length(unique(tl[c(i-1,i)]))==2) breaks<-c(breaks,pbp[i-1])
  return(breaks)
 }
 
 basename<-x$basename
 ll <- readLines(paste(basename,".post",sep=""))
 breaks<-NULL
 for(i in 1:length(ll)) breaks<-c(breaks,makebplist(ll,i))
 cat(paste("The 95% credible interval for a single break point between nucleotide number",leftedge,"and nucleotide number",rightedge,"is:\n"))
 return(quantile(breaks[breaks>leftedge&breaks<rightedge],probs=c(.025,.975)))

}
