sim.codon <-
function(n.taxa,n.codons,reps,tree,template.dat,datfile,freqs=c(0.6,0.3,0.1),omega=c(0.1,0.8,3.2),kap=2){
   tmp<-NA
   while(length(tmp)<(n.taxa+6)){	# added this in for when it screws up
      make.dat(n.taxa,n.codons,reps,tree,template.dat,datfile,freqs,omega,kap)
      command<-paste("/homes/peterchi/bin/evolverNSsites 6",datfile)
      system(command,ignore.stdout=T)
      tmp<-readLines("mc.paml")		# Not really happy about this
      if(length(tmp)<(n.taxa+6)){
         print("would-be error")
      }
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
   return(all.aligns)
}
