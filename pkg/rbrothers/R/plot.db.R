"plot.db" <- function(x,makepic="FALSE",ext="png",seetrees="FALSE",numplot=4,type="unrooted",linewidth=1,colors=NULL,threshold=0) {
if(ext!="png"&ext!="pdf"){
 print("ext must be pdf or png")
 return()
}
library(ape)

basename<-x[[1]]
probs<-x[[2]]
profile<-x[[3]]
treenum<-length(x[[6]])
trees<-x[[6]]

miny<-max(treenum,2)
if(!is.null(colors)){palette(colors)}

columns<-c(2:(dim(probs)[2]))
if(threshold>0){
 columns<-NULL
 for(i in 2:(dim(probs)[2])){
  if((sum(probs[,i]>threshold))>0) columns<-c(columns,i)
 }
 numplot=length(columns)
 miny<-numplot
}

#### plot some probabilities

#### first plots written to file
if(makepic=="TRUE"){
if(ext=="png"){png(paste(basename,".plot1.",ext,sep=""))}
if(ext=="pdf"){pdf(paste(basename,".plot1.",ext,sep=""))}
if(seetrees=="FALSE"){par(mfrow=c(2,1))}
if(seetrees=="TRUE"){
 
 miny2<-max(min(treenum,numplot),1)
 rows<-ceiling(sqrt(miny2))
 cols<-ceiling(miny2/(ceiling(sqrt(miny2))))
 if(miny2==2){
  rows<-1
  cols<-2
 }
 if(threshold>0) miny2<-miny
 layout(matrix(c(rep(1,cols),seq(1:(rows*cols))+1), rows+1, cols, byrow = TRUE))

}

plot(probs$V1,probs[,columns[1]],ylab="Tree Topology Probability",xlab="Nucleotide Position",main=paste(basename,"Recombination Analysis"),type="l",col=1,ylim=c(0,1),lwd=linewidth)
if(numplot>1){for(i in 2:miny){
 lines(probs$V1,probs[,columns[i]],type="l",col=i,lwd=linewidth)}
}
if(numplot==1) lines(probs$V1,probs[,columns[1]],type="l",col=1,lwd=linewidth)
if(seetrees=="FALSE"){plot(profile[,1],profile[,10],type="l",xlab="Nucleotide Position",ylab="Break-Point Probability")}
if(seetrees=="TRUE"){for(i in 1:miny2){plot.phylo(trees[[columns[i]-1]],edge.width=3,edge.col=i,cex=1,type=type)}}
dev.off()
}


if(seetrees=="FALSE"){par(mfrow=c(2,1))}
if(seetrees=="TRUE"){
 miny2<-max(min(treenum,numplot),1)
 rows<-ceiling(sqrt(miny2))
 cols<-ceiling(miny2/(ceiling(sqrt(miny2))))
 if(miny2==2){
  rows<-1
  cols<-2
 }
 if(threshold>0) miny2<-miny
 layout(matrix(c(rep(1,cols),seq(1:(rows*cols))+1), rows+1, cols, byrow = TRUE))
}
plot(probs$V1,probs[,columns[1]],ylab="Tree Topology Probability",xlab="Nucleotide Position",main=paste(basename,"Recombination Analysis"),type="l",col=1,ylim=c(0,1),lwd=linewidth)
if(numplot>1){
 for(i in 2:miny){ 
  lines(probs$V1,probs[,columns[i]],type="l",col=i,lwd=linewidth)}
}
if(numplot==1) lines(probs$V1,probs[,columns[1]],type="l",col=1,lwd=linewidth)
if(seetrees=="FALSE"){plot(profile[,1],profile[,10],type="l",xlab="Nucleotide Position",ylab="Break-Point Probability")}
if(seetrees=="TRUE"){for(i in 1:miny2){plot.phylo(trees[[columns[i]-1]],edge.width=3,edge.col=i,cex=1,type=type)}}

palette("default")

}

