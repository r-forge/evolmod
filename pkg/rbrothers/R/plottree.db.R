"plottree.db" <- function(x,numplot=12,makepng="FALSE") {

library(ape)

basename<-x[[1]]
probs<-x[[2]]
profile<-x[[3]]
treenum<-length(x[[6]])
trees<-x[[6]]

miny<-max(min(treenum,numplot),2)
palette(rainbow(max(treenum,2)))


if(makepng=="TRUE"){
if(treenum==1){miny<-1}
if(numplot==1){miny<-1}
png(paste(basename,".plot2",".png",sep=""))
par(mfrow=c(ceiling(sqrt(miny)),ceiling(miny/(ceiling(sqrt(miny))))))
for(i in 1:miny){
 plot.phylo(trees[[i]],edge.width=3,edge.col=i)
}
dev.off()
}

par(mfrow=c(ceiling(sqrt(miny)),ceiling(miny/(ceiling(sqrt(miny))))))
for(i in 1:miny){
 plot.phylo(trees[[i]],edge.width=3,edge.col=i)
}


}