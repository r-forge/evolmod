calc.null.pboot <-
function(data,B,l,m,exhaustive=F,codeml="codeml",evolver="evolverNSsites"){

   IDnum<-which(!is.na(suppressWarnings(as.numeric(substring(data$nam,1,1)))))
   if(length(IDnum)>0){
      for(i in 1:length(IDnum)){
         data$nam[IDnum[i]]<-paste("N",data$nam[IDnum[i]],sep="")
      }
   }		# append "N" in front of taxa names that start with a number; codeml tanks otherwise

   distance<-dist.dna(as.DNAbin(data),model="F84")*3
   ols.tree<-nnls.tree(distance,bionj(distance),trace=0)
   params<-paml.codeml(data,ols.tree,codeml)
   data<-as.matrix.alignment(data)
   length.codons<-dim(data)[2]/3
   n.taxa<-dim(data)[1]
   Dss.null<-rep(NA,B)


   for(i in 1:B){
      align.taxa<-sim.codon(n.taxa,length.codons,1,ols.tree,freqs=params$freqs,omega=params$omega,kap=params$kap,evolver)[[1]]  
      data<-t(matrix(as.vector(unlist(apply(align.taxa,1,strsplit,""))),ncol=nrow(align.taxa))) # for calc.dss
      data.tmp<-ape::as.alignment(data)

      Dss.null[i]<-max(calc.Dss(data.tmp,l,m,exhaustive))
   }
   return(Dss.null[order(Dss.null)])
}


