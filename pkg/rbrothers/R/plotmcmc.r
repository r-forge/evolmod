plotmcmc<-function(x){
library(coda)
par(mfrow=c(2,1))
traceplot(mcmc(x$MCMC_log_likelihood),ylab="log likelihood",main="MCMC trace plot of log likelihood")
autocorr.plot(x$MCMC_log_likelihood,auto.layout=FALSE)
}
