make.dat <-
function(n.seq,n.codons,reps,tree,freqs,omega,kap){
   ### Create temporary directory
   current.dir<-getwd()
   temp.dir<-paste(current.dir,paste("/paml.evolver.",as.character(round(runif(1,0,1000000000))),sep=""),sep="")
   while(file.exists(temp.dir)){
      temp.dir<-paste(current.dir,paste("/paml.evolver.",as.character(round(runif(1,0,1000000000))),sep=""),sep="")   
   }
   dir.create(temp.dir)
   setwd(temp.dir)

   template<-readLines(system.file("paml/template.dat",package="synDss"))
   template[2]<-round(runif(1,1,1000000))
   template[3]<-paste(n.seq,n.codons,reps,sep=" ")
   template[6]<-tree
   if(length(freqs)!=length(omega)){
      stop("Number of elements in freqs and omega parameters should be equal")
   }
   template[8]<-length(freqs)
   template[9]<-paste("  ",freqs[1],"     ",freqs[2],"    ",freqs[3],sep="")
   template[10]<-paste("  ",omega[1],"     ",omega[2],"    ",omega[3],sep="")
   template[13]<-kap
   writeLines(template,"paml.dat")
   return(temp.dir)
}



