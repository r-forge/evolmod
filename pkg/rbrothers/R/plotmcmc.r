plotmcmc<-function(x,mode="trace"){
library(coda)
if(mode!="trace"&mode!="auto"){
print("mode must be trace or auto")
return()
}
par(mfrow=c(1,1))
if(mode=="trace"){traceplot(mcmc(x$MCMC_log_likelihood),ylab="log likelihood",main="MCMC trace plot of log likelihood")}
if(mode=="auto"){autocorr.plot(x$MCMC_log_likelihood,auto.layout=FALSE)}
}
