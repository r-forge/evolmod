
"recombinationplots"<-function(x,visualize="FALSE"){
basename<-x$basename
library(ape)
probs<-read.table(paste(basename,".topoprob",sep=""))
profile<-read.table(paste(basename,".profile",sep=""))

#### creates a vector of all the tip names culled from the .phy file
P<-as.integer(scan(paste(basename,".phy",sep=""),what='raw',nlines=1)[1])
names<-NULL
for(i in 1:P){
names[i]<-(scan(paste(basename,".phy",sep=""),what='raw',nlines=1,skip=i)[1])
}

#### read in the trees
treenum<-dim(probs)[2]-1
trees<-read.tree(paste(basename,".tree",sep=""),tree.names=seq(1:treenum),keep.multi=TRUE)

#### relabel the tips with the names from the .phy file
for(i in 1:treenum){
 for(j in 1:P){
  trees[[i]]$tip.label[j]<-names[as.integer(trees[[i]]$tip.label[j])+1]
 }
}

miny<-max(treenum,2)
palette(rainbow(miny))

#### plot some probabilities

#### first plots written to file

png(paste(basename,".plot1",".png",sep=""))
par(mfrow=c(2,1))
plot(probs$V1,probs$V2,ylab="Tree Topology Probability",xlab="Nucleotide Position",main=paste(basename,"Recombination Analysis"),type="l",col=1,ylim=c(0,1))
if(treenum>1){for(i in 3:(miny+1)){lines(probs$V1,probs[,i],type="l",col=(i-1))}}
plot(profile[,1],profile[,10],type="l",xlab="Nucleotide Position",ylab="Break-Point Probability")
dev.off()
#print(paste(basename,".plot1",".png ","created (hopefully).",sep=""))

if(visualize=="TRUE"){
par(mfrow=c(2,1))
plot(probs$V1,probs$V2,ylab="Tree Topology Probability",xlab="Nucleotide Position",main=paste(basename,"Recombination Analysis"),type="l",col=1,ylim=c(0,1))
if(treenum>1){for(i in 3:(miny+1)){lines(probs$V1,probs[,i],type="l",col=(i-1))}}
plot(profile[,1],profile[,10],type="l",xlab="Nucleotide Position",ylab="Break-Point Probability")
}
}

