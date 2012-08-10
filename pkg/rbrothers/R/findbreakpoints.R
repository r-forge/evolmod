"findbreakpoints"<-function(basename){
 probs<-read.table(paste(basename,".topoprob",sep=""))
 toptree<-findmax(probs,1)
 breaks<-NULL
 k<-0
 for(i in 2:dim(probs)[1]){
  nexty<-findmax(probs,i)
  if(nexty!=toptree){
   toptree<-nexty
   k<-k+1
   breaks[k]<-i
  }
 }
 write(breaks, file = paste(basename,".break",sep=""))
 return(breaks)
}
