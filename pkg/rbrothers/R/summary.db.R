"summary.db"<-function(x, ..){
cat(paste("DualBrothers output for",x[[1]]),"\n")
cat(paste(x[[4]],"sequences of length",dim(x[[2]])[1]),"\n")
if(length(x[[5]])>1){cat(length(x[[5]]),"breakpoints:",c(x[[5]]),"\n")}
if(length(x[[5]])==1){cat(length(x[[5]]),"breakpoint:",c(x[[5]]),"\n")}
if(length(x[[5]])==0){cat("no breakpoints","\n")}
if(length(x[[6]])>1){cat(paste(length(x[[6]]),"trees"),"\n")}
if(length(x[[6]])==1){cat(paste(length(x[[6]]),"tree"),"\n")}
}