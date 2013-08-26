calc.null.pboot.syn <-
function(data,B,l,m,syn.matrix,exhaustive=F){
   length.codons<-dim(data)[2]/3
   n.taxa<-dim(data)[1]
   Dss.null<-rep(NA,B)
   distance<-dist.dna(as.DNAbin(as.alignment(data)),model="F84")*3
   ols.tree<-optim.phylo.ls.all(distance)
   id<-.Random.seed[1]
   ctl.file<-paste(paste("pboot",id,sep=""),".ctl",sep="")
   params<-codeml.M3(id,data,ctl.file,ols.tree)   
   for(i in 1:B){
      align.taxa<-sim.codon(n.taxa,length.codons,1,write.tree(ols.tree),"template.dat","boot.dat",freqs=params$freqs,omega=params$omega,kap=params$kap)[[1]]  
      data.tmp<-read.phylosim(align.taxa)

      Dss.null[i]<-max(calc.Dss.syn(data.tmp,l,m,syn.matrix,exhaustive))
   }
   return(Dss.null[order(Dss.null)])
}
