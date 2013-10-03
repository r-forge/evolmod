"readdb"<-function(basename){
probs<-read.table(paste(basename,".topoprob",sep=""))
profile<-read.table(paste(basename,".profile",sep=""))
P<-as.integer(scan(paste(basename,".phy",sep=""),what='raw',nlines=1)[1])
names<-NULL
for(i in 1:P){
names[i]<-(scan(paste(basename,".phy",sep=""),what='raw',nlines=1,skip=i)[1])
}
breaks<-findbreakpoints(basename)
treenum<-dim(probs)[2]-1
trees<-read.tree(paste(basename,".tree",sep=""),tree.names=seq(1:treenum),keep.multi=TRUE)
for(i in 1:treenum){
 for(j in 1:P){
  trees[[i]]$tip.label[j]<-names[as.integer(trees[[i]]$tip.label[j])+1]
 }
}

#### read the command file
cmdf<-read.table(paste(basename,".cmdfile",sep=""))
####


#cmd<-scan(paste(basename,".cmdfile",sep=""),what='raw',nlines=32)
#for(i in 1:(length(cmd)/2)){
# if(cmd[2*(i-1)+1]=="par_lambda:"){par_lambda<-as.numeric(cmd[2*i])}
# if(cmd[2*(i-1)+1]=="top_lambda:"){top_lambda<-as.numeric(cmd[2*i])}
# if(cmd[2*(i-1)+1]=="sigma_alpha:"){sigma_alpha<-as.numeric(cmd[2*i])}
# if(cmd[2*(i-1)+1]=="sigma_mu:"){sigma_mu<-as.numeric(cmd[2*i])}
# if(cmd[2*(i-1)+1]=="top_breaks:"){top_breaks<-as.numeric(cmd[2*i])}
# if(cmd[2*(i-1)+1]=="par_breaks:"){par_breaks<-as.numeric(cmd[2*i])}
#}


####
#### read the .post file to find the number of topology breakpoints and log likelihoods
####

lineiloglik<-function(linei){return(sum(as.numeric(linei[5+8*(0:(as.numeric(linei[2])-1))])))}

lineinumbrkpnts<-function(linei){
 numchgpnts<-as.numeric(linei[2])
 ts1<-linei[4]
 returnme<-0
 if(numchgpnts>1) for(i in 2:numchgpnts) if(!linei[4+8*(i-2)]==linei[4+8*(i-1)]){returnme<-returnme+1}
 return(returnme)
}

loglikes<-NULL
numbrkpnts<-NULL
ll <- readLines(paste(basename,".post",sep=""))

for(i in 1:length(ll)){
 ii = unlist(strsplit(ll[[i]]," "))
 ii = ii[which(sapply(ii,nchar,USE.NAMES=FALSE)!=0)]
 loglikes[i]<-lineiloglik(ii)
 numbrkpnts[i]<-lineinumbrkpnts(ii)
}

####
####
####

tbr<-list(basename,probs,profile,P,breaks,trees,cmdf,loglikes,numbrkpnts)
class(tbr)<-"db"
names(tbr)<-c("basename","TopologyProfile","EPProfile","numberofsequences","breakpoints","trees","command_file","MCMC_log_likelihood","MCMC_number_of_break_points")
return(tbr)
}