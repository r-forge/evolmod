"summary.db"<-function(x, ..){
cat(paste("DualBrothers output for",x[[1]]),"\n")
cat(paste(x[[4]],"sequences of length",dim(x[[2]])[1]),"\n")
#if(length(x[[5]])>1){cat(length(x[[5]]),"breakpoints:",c(x[[5]]),"\n")}
#if(length(x[[5]])==1){cat(length(x[[5]]),"breakpoint:",c(x[[5]]),"\n")}
#if(length(x[[5]])==0){cat("no breakpoints","\n")}
cat(paste("average number of breakpoints in the posterior:",round(mean(x$MCMC_number_of_break_points),1)),"\n")
if(length(x[[6]])>1){cat(paste(length(x[[6]]),"trees"),"\n\n")}
if(length(x[[6]])==1){cat(paste(length(x[[6]]),"tree"),"\n\n")}
cat(paste("Prior parameters:"),"\n")
cat(paste("par_lambda:",x$par_lambda),"\n")
cat(paste("top_lambda:",x$top_lambda),"\n")
cat(paste("sigma_alpha:",x$sigma_alpha),"\n")
cat(paste("sigma_mu:",x$sigma_mu),"\n")
}