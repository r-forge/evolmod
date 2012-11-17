"summary.db"<-function(x, ..){

### parsing command file elements:

for(i in 1:(dim(x$command_file)[1])){
 if(x$command_file[i,1]=="par_lambda:"){par_lambda<-as.numeric(as.character(x$command_file[i,2]))}
 if(x$command_file[i,1]=="top_lambda:"){top_lambda<-as.numeric(as.character(x$command_file[i,2]))}
 if(x$command_file[i,1]=="length:"){length<-as.numeric(as.character(x$command_file[i,2]))}
 if(x$command_file[i,1]=="burnin:"){burnin<-as.numeric(as.character(x$command_file[i,2]))}
 if(x$command_file[i,1]=="subsample:"){subsample<-as.numeric(as.character(x$command_file[i,2]))}
}

### Bayes factor calculation:
#bf<-(sum(x$MCMC_number_of_break_points>0)/sum(x$MCMC_number_of_break_points==0))/(exp(top_lambda)-1)
#if(bf!=Inf) bf<-signif(bf,2)
###

### a slightly different Bayes factor calculation:
Y<-1000*(exp(top_lambda)-1)
rm<-Y/(1+Y)
###

cat(paste("DualBrothers output for",x[[1]]),"\n")
cat(paste(x[[4]],"sequences of length",dim(x[[2]])[1]),"\n\n")

cat(paste("MCMC settings:"),"\n")
cat(paste("length of the MCMC chain:",length),"\n")
cat(paste("burn-in length:",burnin),"\n")
cat(paste("subsample frequency:",subsample),"\n\n")

cat(paste("Prior parameters:"),"\n")
cat(paste("prior mean number of substitution process change points:",par_lambda),"\n")
cat(paste("prior mean number of topology change points:",top_lambda),"\n\n")

cat(paste("average number of breakpoints in the posterior:",round(mean(x$MCMC_number_of_break_points),1)),"\n")
cat(paste("The posterior probability of at least one breakpoint was",round(sum(x$MCMC_number_of_break_points>0)/length(x$MCMC_number_of_break_points),3)," (>",round(rm,3)," required for a Bayes factor > 1000).","\n"))
#if(sum(x$MCMC_number_of_break_points==0)==0){cat(paste("All of the posterior samples contained at least one break point."),"\n")}
#if(bf!=Inf){cat(paste("The Bayes factor in favor of recombination is ",bf),"\n")}
if(length(x[[6]])>1){cat(paste(length(x[[6]]),"trees considered"),"\n\n")}
if(length(x[[6]])==1){cat(paste(length(x[[6]]),"tree considered"),"\n\n")}

}