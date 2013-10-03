"plottree.db" <- function(x,numplot=4,makepic="FALSE",ext="png",type="unrooted",colors=NULL,threshold=0) {

basename<-x[[1]]
probs<-x[[2]]
profile<-x[[3]]
treenum<-length(x[[6]])
trees<-x[[6]]

miny<-max(min(treenum,numplot),1)
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

if(makepic=="TRUE"){
if(treenum==1){miny<-1}
if(numplot==1){miny<-1}
if(ext=="png"){png(paste(basename,"plot2",".png",sep=""))}
if(ext=="pdf"){pdf(paste(basename,"plot2",".pdf",sep=""))}
par(mfrow=c(ceiling(sqrt(miny)),ceiling(miny/(ceiling(sqrt(miny))))))
for(i in 1:miny){
 plot.phylo(trees[[columns[i]-1]],edge.width=3,edge.color=i,type=type)
}
dev.off()
}

par(mfrow=c(ceiling(sqrt(miny)),ceiling(miny/(ceiling(sqrt(miny))))))
for(i in 1:miny){
 plot.phylo(trees[[columns[i]-1]],edge.width=3,edge.color=i,type=type)
}


}