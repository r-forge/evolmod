calc.null.pboot <-
function(data,B,l,m,exhaustive=F){
   distance<-dist.dna(as.DNAbin(data),model="F84")*3
   ols.tree<-optim.phylo.ls.all(distance)
   data<-as.matrix.alignment(data)
   length.codons<-dim(data)[2]/3
   n.taxa<-dim(data)[1]
   Dss.null<-rep(NA,B)

   params<-paml.codeml(data,ols.tree)

   for(i in 1:B){
      align.taxa<-sim.codon(n.taxa,length.codons,1,write.tree(ols.tree),"template.dat","boot.dat",freqs=params$freqs,omega=params$omega,kap=params$kap)[[1]]  
      data<-t(matrix(as.vector(unlist(apply(align.taxa,1,strsplit,""))),ncol=nrow(align.taxa))) # for calc.dss
      data.tmp<-ape::as.alignment(data)

      Dss.null[i]<-max(calc.Dss(data.tmp,l,m,exhaustive))
   }
   return(Dss.null[order(Dss.null)])
}


