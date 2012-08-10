"topologyprofile"<-function(basename){
   seqlength<-scan(paste(basename,".phy",sep=""),nlines=1,quiet= TRUE)[2]
   tps <- .jnew("TopologyProfile", paste(basename,".post",sep=""),paste(basename,".topoprob",sep=""),paste(basename,".tree",sep=""),as.integer(seqlength))
   #system(paste("java TopologyProfile ",basename,".post ",basename,".topoprob ",basename,".tree ",seqlength,sep=""))
}
