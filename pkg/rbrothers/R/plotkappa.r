plotkappa<-function(x){
 basename<-x$basename
 epp<-read.table(paste(basename,".profile",sep=""))
 plot(epp[,2]~epp[,1],type="l",col="red",ylim=c(min(epp[,4]),max(epp[,6])),xlab="position along sequence",ylab="Kappa value",main="Posterior Kappa Values")
 lines(epp[,4]~epp[,1],lty=2,col="green")
 lines(epp[,5]~epp[,1],lty=1,col="blue")
 lines(epp[,6]~epp[,1],lty=2,col="purple")
 legend("topleft",legend=c("mean","median","2.5% quantile","97.5% quantile"),lty=c(1,1,2,2),col=c("red","blue","green","purple"))
}
