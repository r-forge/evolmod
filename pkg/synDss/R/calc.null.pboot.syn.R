calc.null.pboot.syn <-
function(data,B,l,m,syn.matrix,exhaustive=F,codeml="codeml",evolver="evolverNSsites 6"){
   distance<-dist.dna(as.DNAbin(data),model="F84")*3
   ols.tree<-nnls.tree(distance,bionj(distance),trace=0)
   params<-paml.codeml(data,ols.tree,codeml)   
   data<-make.mj.align(as.matrix.alignment(data))
   length.codons<-dim(data)[2]/3
   n.taxa<-dim(data)[1]
   Dss.null<-rep(NA,B)
   for(i in 1:B){
      align.taxa<-sim.codon(n.taxa,length.codons,1,ols.tree,freqs=params$freqs,omega=params$omega,kap=params$kap,evolver)[[1]]  
      data<-t(matrix(as.vector(unlist(apply(align.taxa,1,strsplit,""))),ncol=nrow(align.taxa)))
      data.tmp<-ape::as.alignment(data)

      Dss.null[i]<-max(calc.Dss.syn(data.tmp,l,m,syn.matrix,exhaustive))
   }
   return(Dss.null[order(Dss.null)])
}


