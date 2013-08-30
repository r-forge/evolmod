calc.null.pboot.syn <-
function(data,B,l,m,syn.matrix,exhaustive=F){
   distance<-dist.dna(as.DNAbin(data),model="F84")*3
   ols.tree<-optim.phylo.ls.all(distance)
   params<-paml.codeml(data,ols.tree)   
   data<-as.matrix.alignment(make.mj.align(data))   # 8/29 not yet tested
   length.codons<-dim(data)[2]/3
   n.taxa<-dim(data)[1]
   Dss.null<-rep(NA,B)
   for(i in 1:B){
      align.taxa<-sim.codon(n.taxa,length.codons,1,write.tree(ols.tree),freqs=params$freqs,omega=params$omega,kap=params$kap)[[1]]  
      data.tmp<-make.mj.align(align.taxa)

      Dss.null[i]<-max(calc.Dss.syn(data.tmp,l,m,syn.matrix,exhaustive))
   }
   return(Dss.null[order(Dss.null)])
}


