"plottree.db" <- function(x,numplot=4,makepic="FALSE",ext="png",type="unrooted") {

library(ape)

basename<-x[[1]]
probs<-x[[2]]
profile<-x[[3]]
treenum<-length(x[[6]])
trees<-x[[6]]

miny<-max(min(treenum,numplot),1)
#palette(rainbow(max(treenum,2)))


if(makepic=="TRUE"){
if(treenum==1){miny<-1}
if(numplot==1){miny<-1}
if(ext=="png"){png(paste(basename,".plot2",".png",sep=""))}
if(ext=="pdf"){pdf(paste(basename,".plot2",".pdf",sep=""))}
par(mfrow=c(ceiling(sqrt(miny)),ceiling(miny/(ceiling(sqrt(miny))))))
for(i in 1:miny){
 plot.phylo(trees[[i]],edge.width=3,edge.col=i,type=type)
}
dev.off()
}

par(mfrow=c(ceiling(sqrt(miny)),ceiling(miny/(ceiling(sqrt(miny))))))
for(i in 1:miny){
 plot.phylo(trees[[i]],edge.width=3,edge.col=i,type=type)
}


}