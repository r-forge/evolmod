"plot.db" <- function(x,makepic="FALSE",ext="png",seetrees="FALSE",numplot=4,type="unrooted") {
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
#palette(rainbow(miny))

#### plot some probabilities

#### first plots written to file

if(makepic=="TRUE"){
if(ext=="png"){png(paste(basename,".plot1.",ext,sep=""))}
if(ext=="pdf"){pdf(paste(basename,".plot1.",ext,sep=""))}
if(seetrees=="FALSE"){par(mfrow=c(2,1))}
if(seetrees=="TRUE"){
 miny2<-max(min(treenum,numplot),1)
 rows<-ceiling(sqrt(miny2))
 columns<-ceiling(miny2/(ceiling(sqrt(miny2))))
 if(miny2==2){
  rows<-1
  columns<-2
 }
 layout(matrix(c(rep(1,columns),seq(1:(rows*columns))+1), rows+1, columns, byrow = TRUE))
}
plot(probs$V1,probs$V2,ylab="Tree Topology Probability",xlab="Nucleotide Position",main=paste(basename,"Recombination Analysis"),type="l",col=1,ylim=c(0,1))
if(treenum>1){for(i in 3:(miny+1)){lines(probs$V1,probs[,i],type="l",col=(i-1))}}
if(seetrees=="FALSE"){plot(profile[,1],profile[,10],type="l",xlab="Nucleotide Position",ylab="Break-Point Probability")}
if(seetrees=="TRUE"){for(i in 1:miny2){plot.phylo(trees[[i]],edge.width=3,edge.col=i,cex=1,type=type)}}
dev.off()
}

miny<-max(treenum,2)
if(seetrees=="FALSE"){par(mfrow=c(2,1))}
if(seetrees=="TRUE"){
 miny2<-max(min(treenum,numplot),1)
 rows<-ceiling(sqrt(miny2))
 columns<-ceiling(miny2/(ceiling(sqrt(miny2))))
 if(miny2==2){
  rows<-1
  columns<-2
 }
 layout(matrix(c(rep(1,columns),seq(1:(rows*columns))+1), rows+1, columns, byrow = TRUE))
}
plot(probs$V1,probs$V2,ylab="Tree Topology Probability",xlab="Nucleotide Position",main=paste(basename,"Recombination Analysis"),type="l",col=1,ylim=c(0,1))
if(treenum>1){for(i in 3:(miny+1)){lines(probs$V1,probs[,i],type="l",col=(i-1))}}
if(seetrees=="FALSE"){plot(profile[,1],profile[,10],type="l",xlab="Nucleotide Position",ylab="Break-Point Probability")}
if(seetrees=="TRUE"){for(i in 1:miny2){plot.phylo(trees[[i]],edge.width=3,edge.col=i,cex=1,type=type)}}

}