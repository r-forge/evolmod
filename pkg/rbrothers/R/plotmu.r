plotmu<-function(x){
basename<-x$basename
epp<-read.table(paste(basename,".profile",sep=""))
plot(epp[,3]~epp[,1],type="l",col="red",ylim=c(min(epp[,7]),max(epp[,9])),xlab="position along sequence",ylab="Mu value",main="Posterior Mu Values")
lines(epp[,7]~epp[,1],lty=2,col="green")
lines(epp[,8]~epp[,1],lty=1,col="blue")
lines(epp[,9]~epp[,1],lty=2,col="purple")
legend("topleft",legend=c("mean","median","2.5% quantile","97.5% quantile"),lty=c(1,1,2,2),col=c("red","blue","green","purple"))
}
