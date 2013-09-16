paml.codeml <-
function(data, tree, call="codeml"){
   
  ### Create temporary directory and files.
   current.dir<-getwd()
   temp.dir<-paste(current.dir,paste("/paml.codeml.",as.character(round(runif(1,0,10000000))),sep=""),sep="")
   if(file.exists(temp.dir)){
      unlink(temp.dir, recursive = TRUE)
   } else{
      dir.create(temp.dir)
   }


   temp.file.tree <- paste(temp.dir, "/codeml.trees", sep = "")
   temp.file.nuc <- paste(temp.dir, "/codeml.nuc", sep = "")
   temp.file.control <- paste(temp.dir, "/template.ctl", sep = "")
   baseml.file.names <- c("codeml.trees", "codeml.nuc", "codeml.ctl")

   temp.file.stdout <- paste(temp.dir, "/stdout", sep = "")

   # write sequence and tree files
   write.dna(as.matrix.alignment(data),temp.file.nuc,format="fasta")
   write.tree(tree,temp.file.tree)

   # read in template control file
   template<-readLines(system.file("paml/template.ctl", package="synDss"))

   # modify template file for current project
   template[1]<-"      seqfile = codeml.nuc"
   template[2]<-"     treefile = codeml.trees"
   template[4]<-"      outfile = stdout"

   # write new control file
   writeLines(template,temp.file.control)

   # run codeml
   setwd(temp.dir)
   command<-paste(call,"template.ctl",sep=" ")
   system(command,ignore.stdout=T,show.output.on.console=F)
   setwd(current.dir)

   # read and parse output
   out<-readLines(temp.file.stdout)
   kap<-as.numeric(unlist(strsplit(out[grep("kap",out)],"[[:blank:]]+"))[4])
   freqs<-as.numeric(unlist(strsplit(out[grep("p:   ",out)],"[[:blank:]]+"))[c(2,3,4)])
   omegas<-as.numeric(unlist(strsplit(out[grep("w:   ",out)],"[[:blank:]]+"))[c(2,3,4)])
   output<-list(kap,freqs,omegas)
   names(output)<-c("kap","freqs","omegas")

   # cleanup temporary files
   unlink(temp.dir, recursive=TRUE)

   return(output)

}



