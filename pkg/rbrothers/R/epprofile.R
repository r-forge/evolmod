"epprofile"<-function(basename){
  seqlength<-scan(paste(basename,".phy",sep=""),nlines=1,quiet= TRUE)[2]
  out <- .jnew("EPProfile", paste(basename,".post",sep=""),paste(basename,".profile",sep=""),as.integer(seqlength))
  #system(paste("java EPProfile ",basename,".post ",basename,".profile ",seqlength,sep=""))
}
