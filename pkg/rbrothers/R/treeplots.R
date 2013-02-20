
"treeplots"<-function(x,numplot=12,visualize="FALSE"){
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

#### maximum number of trees to plot can be changed here (default 12)
miny<-max(min(treenum,numplot),1)
palette(rainbow(max(treenum,2)))
#palette(rainbow(miny))

#### plot the corresponding trees

#### second plots written to file
if(treenum==1){miny<-1}
if(numplot==1){miny<-1}
png(paste(basename,".plot2",".png",sep=""))
par(mfrow=c(ceiling(sqrt(miny)),ceiling(miny/(ceiling(sqrt(miny))))))
for(i in 1:miny){
 plot.phylo(trees[[i]],edge.width=3,edge.col=i)
}
dev.off()

if(visualize=="TRUE"){
par(mfrow=c(ceiling(sqrt(miny)),ceiling(miny/(ceiling(sqrt(miny))))))
for(i in 1:miny){
 plot.phylo(trees[[i]],edge.width=3,edge.col=i)
}
}
}
