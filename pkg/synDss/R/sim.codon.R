sim.codon <-
function(n.taxa,n.codons,reps,tree,freqs=c(0.6,0.3,0.1),omega=c(0.1,0.8,3.2),kap=2,call="evolverNSsites"){
   tree<-write.tree(tree)
   # need to add error checking, or else it will create directories and crash before it can remove them
   current.dir<-getwd()
   tmp<-NA
   temp.dir<-make.dat(n.taxa,n.codons,reps,tree,freqs,omega,kap) # write paml.dat and return the dir it is in
   setwd(temp.dir)	
   command<-paste(call,"6 paml.dat")
   err1<-try(system(command,ignore.stdout=T,show.output.on.console=F))
   if(err1!=0){
      setwd(current.dir)
      unlink(temp.dir, recursive=TRUE)
      stop("error in evolver call")
   }
   err2<-is.element("mc.paml",list.files())
   if(!err2){
      setwd(current.dir)
      unlink(temp.dir, recursive=TRUE)
      stop("error in reading mc.paml file")
   } else {
      tmp<-readLines("mc.paml")
   }
   err3<-length(tmp)<(n.taxa+6)
   if(err3){
      setwd(current.dir)
      unlink(temp.dir, recursive=TRUE)
      stop("error in mc.paml file")
   }

   starts<-intersect(grep(n.taxa,tmp),grep(n.codons*3,tmp))+1	# the rows that have both n.taxa and n.codons
   all.aligns<-vector("list",reps)
   for(i in 1:reps){
      data.temp<-matrix(NA,ncol=n.codons,nrow=n.taxa)
      rownames.temp<-rep(NA,n.taxa)
      for(j in 1:n.taxa){
         temp<-strsplit(tmp[(starts[i]+j)]," ")[[1]]
         temp2<-temp[!(temp=="")]
         data.temp[j,]<-temp2[-1]
         rownames.temp[j]<-as.character(temp2[1])
      }
   rownames(data.temp)<-rownames.temp
   colnames(data.temp)<-seq(1:dim(data.temp)[2])
   all.aligns[[i]]<-data.temp
   }

   setwd(current.dir)

   # cleanup temporary files
   unlink(temp.dir, recursive=TRUE)

   return(all.aligns)
}



