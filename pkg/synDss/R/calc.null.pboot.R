calc.null.pboot <-
function(data,B,l,m,exhaustive=F){
   length.codons<-dim(data)[2]/3
   n.taxa<-dim(data)[1]
   Dss.null<-rep(NA,B)
   distance<-dist.dna(as.DNAbin(as.alignment(data)),model="F84")*3
   ols.tree<-optim.phylo.ls.all(distance)
   id<-abs(.Random.seed[3])%%10^3	# just get some 3-digit identification of the seed
   ctl.file<-paste(paste("pboot",id,sep=""),".ctl",sep="")
   params<-codeml.M3(id,data,ctl.file,ols.tree)

   for(i in 1:B){
      align.taxa<-sim.codon(n.taxa,length.codons,1,write.tree(ols.tree),"template.dat","boot.dat",freqs=params$freqs,omega=params$omega,kap=params$kap)[[1]]  
      data<-t(matrix(as.vector(unlist(apply(align.taxa,1,strsplit,""))),ncol=nrow(align.taxa))) # for calc.dss
      data.tmp<-as.alignment(data)

      Dss.null[i]<-max(calc.Dss.all(data.tmp,l,m,exhaustive))
   }
   return(Dss.null[order(Dss.null)])
}


