# in progress 8/14 -- also need to do calc.Dss.syn and bootstraps.
calc.Dss <-
function(data,l,m,exhaustive=FALSE){
   data<-as.DNAbin(data)
   length<-dim(data)[2]
   n.seq<-dim(data)[1]
   if(exhaustive & n.seq>6){
      stop("Exhaustive search not allowed with greater than 6 sequences")
   }
   n.windows<-((length-(2*l))/m+1)
   d.bar<-mean(dist.dna(data,model="JC69"))
   Dss.all<-rep(NA,n.windows)
   for(i in 1:n.windows){
      start.a<-1+m*(i-1)
      end.a<-m*(i-1)+l
      start.b<-1+l+m*(i-1)
      end.b<-2*l+m*(i-1)
      
      dist1<-dist.dna(data[,start.a:end.a],model="JC69")
      w1<-mean(dist1)
      dist1.stan<-dist1*(d.bar/w1)

      dist2<-dist.dna(data[,start.b:end.b],model="JC69")
      w2<-mean(dist2)
      dist2.stan<-dist2*(d.bar/w2)

      if(w1!=0 & w2!=0){		# skip it if either window has identical sequences...
         T1i<-ifelse(exhaustive,optim.phylo.ls.all(dist1.stan),nnls.tree(dist1.stan,bionj(dist1.stan),trace=0))

         SSa.F<-attr(T1i,"RSS")
         SSb.F<-attr(nnls.tree(as.matrix(dist2.stan),T1i,trace=0),"RSS")

         T2i<-ifelse(exhaustive,optim.phylo.ls.all(dist2.stan),nnls.tree(dist2.stan,bionj(dist2.stan),trace=0))
         SSa.B<-attr(T2i,"RSS")
         SSb.B<-attr(nnls.tree(as.matrix(dist1.stan),T2i,trace=0),"RSS")

         Dss.all[i]<-max((SSb.F-SSa.F),(SSb.B-SSa.B))
      }
      else{
         Dss.all[i]<-0			#... and then just make Dss equal 0
      }
   }
   return(Dss.all)
}
