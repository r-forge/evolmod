"plot.db" <- function(x,makepng="FALSE") {

library(ape)

basename<-x[[1]]
probs<-x[[2]]
profile<-x[[3]]
treenum<-length(x[[6]])

miny<-max(treenum,2)
palette(rainbow(miny))

#### plot some probabilities

#### first plots written to file

if(makepng=="TRUE"){
png(paste(basename,".plot1",".png",sep=""))
par(mfrow=c(2,1))
plot(probs$V1,probs$V2,ylab="Tree Topology Probability",xlab="Nucleotide Position",main=paste(basename,"Recombination Analysis"),type="l",col=1,ylim=c(0,1))
if(treenum>1){for(i in 3:(miny+1)){lines(probs$V1,probs[,i],type="l",col=(i-1))}}
plot(profile[,1],profile[,10],type="l",xlab="Nucleotide Position",ylab="Break-Point Probability")
dev.off()
}


par(mfrow=c(2,1))
plot(probs$V1,probs$V2,ylab="Tree Topology Probability",xlab="Nucleotide Position",main=paste(basename,"Recombination Analysis"),type="l",col=1,ylim=c(0,1))
if(treenum>1){for(i in 3:(miny+1)){lines(probs$V1,probs[,i],type="l",col=(i-1))}}
plot(profile[,1],profile[,10],type="l",xlab="Nucleotide Position",ylab="Break-Point Probability")


}