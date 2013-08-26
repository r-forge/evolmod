make.dat <-
function(n.seq,n.codons,reps,tree,template.dat,file,freqs,omega,kap){
   template<-readLines(template.dat)
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
   writeLines(template,file)
}
