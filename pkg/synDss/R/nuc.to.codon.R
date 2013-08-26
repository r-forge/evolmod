nuc.to.codon <-
function(alignment){
   n.nuc<-dim(alignment)[2]
   n.taxa<-dim(alignment)[1]
   if(length(alignment)/3 == round(length(alignment)/3)){
      codon.align<-matrix(NA,ncol=n.nuc/3,nrow=n.taxa)
      for(i in 1:dim(codon.align)[1]){
         for(j in 1:dim(codon.align)[2]){
            codon.align[i,j]<-paste(alignment[i,(3*j-2)],alignment[i,(3*j-1)],alignment[i,(3*j)],sep="")
         }
      }
   }
   rownames(codon.align)<-rownames(alignment)
   colnames(codon.align)<-seq(1:dim(codon.align)[2])
   return(codon.align)
}
