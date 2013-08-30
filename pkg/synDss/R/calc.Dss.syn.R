calc.Dss.syn <-
function(data,l,m,syn.matrix,exhaustive=FALSE){
   data<-as.matrix.alignment(make.mj.align(data))   # 8/29 not yet tested
   length<-dim(data)[2]
   n.seq<-dim(data)[1]
   if(exhaustive & n.seq>6){
      stop("Exhaustive search not allowed with greater than 6 sequences")
   }
   n.windows<-((length-(2*l))/m+1)
   d.bar<-mean(as.dist(robust.codon.dist.f84(data,syn.matrix,partition=1))) 
   Dss.all<-rep(NA,n.windows)
   for(i in 1:n.windows){
      start.a<-1+m*(i-1)
      end.a<-m*(i-1)+l
      start.b<-1+l+m*(i-1)
      end.b<-2*l+m*(i-1)

      temp1<-try(robust.codon.dist.f84(data[,start.a:end.a],syn.matrix,partition=1),T)
      if(class(temp1)=="try-error"){
         dist1<-as.dist(matrix(0,ncol=dim(data)[1],nrow=dim(data)[1]))
      } else{
         dist1<-as.dist(temp1)
      }
      w1<-mean(dist1,na.rm=T)
      dist1.stan<-dist1*(d.bar/w1)

      temp2<-try(robust.codon.dist.f84(data[,start.b:end.b],syn.matrix,partition=1),T)
      if(class(temp2)=="try-error"){
         dist2<-as.dist(matrix(0,ncol=dim(data)[1],nrow=dim(data)[1]))
      } else{
         dist2<-as.dist(temp2)
      }
      w2<-mean(dist2,na.rm=T)
      dist2.stan<-dist2*(d.bar/w2)

      if(w1!=0 & w2!=0){

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
