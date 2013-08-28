codeml.M3 <-
function(id,data,ctl.name,tree){
   title1<-paste(paste("temp",id,sep=""),".nuc",sep="")
   write.dna(data,title1,format="fasta")

   title2<-paste(paste("temp",id,sep=""),".trees",sep="")
   write.tree(tree,title2)

   title3<-paste(paste("mlc",id,sep=""))

   template<-readLines(system.file("paml/template.ctl", package="synDss"))
   template[1]<-paste("      seqfile = ",title1,sep="")
   template[2]<-paste("     treefile = ",title2,sep="")

   template[4]<-paste("      outfile = ",title3,sep="")

   writeLines(template,ctl.name)
   command<-paste("/homes/peterchi/bin/codeml ",ctl.name,sep="")
   system(command,ignore.stdout=T)

#writeLines(paste("Here is the command:",command,sep=" "),title1)

   out<-readLines(title3)
   kap<-as.numeric(unlist(strsplit(out[grep("kap",out)],"[[:blank:]]+"))[4])
   freqs<-as.numeric(unlist(strsplit(out[grep("p:   ",out)],"[[:blank:]]+"))[c(2,3,4)])
   omegas<-as.numeric(unlist(strsplit(out[grep("w:   ",out)],"[[:blank:]]+"))[c(2,3,4)])
   output<-list(kap,freqs,omegas)
   names(output)<-c("kap","freqs","omegas")

   return(output)
}
