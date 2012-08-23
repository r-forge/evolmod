"readdb"<-function(basename){
library(ape)
probs<-read.table(paste(basename,".topoprob",sep=""))
profile<-read.table(paste(basename,".profile",sep=""))
P<-as.integer(scan(paste(basename,".phy",sep=""),what='raw',nlines=1)[1])
names<-NULL
for(i in 1:P){
names[i]<-(scan(paste(basename,".phy",sep=""),what='raw',nlines=1,skip=i)[1])
}
breaks<-findbreakpoints(basename)
treenum<-dim(probs)[2]-1
trees<-read.tree(paste(basename,".tree",sep=""),tree.names=seq(1:treenum),keep.multi=TRUE)
for(i in 1:treenum){
 for(j in 1:P){
  trees[[i]]$tip.label[j]<-names[as.integer(trees[[i]]$tip.label[j])+1]
 }
}
tbr<-list(basename,probs,profile,P,breaks,trees)
class(tbr)<-"db"
names(tbr)<-c("basename","TopologyProfile","EPProfile","numberofsequences","breakpoints","trees")
return(tbr)
}