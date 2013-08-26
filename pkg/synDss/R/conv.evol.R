conv.evol <-
function(alignment,perc,taxa1,taxa2,map){
   stop.codons<-c("TAA","TAG","TGA")
   nucl<-c("A","G","C","T")
   align.taxa.conv<-alignment

   taxa1.aa<-codon.map[match(alignment[taxa1,],map[,1]),2]
   taxa2.aa<-codon.map[match(alignment[taxa2,],map[,1]),2]

   ind<-which(taxa1.aa!=taxa2.aa)
   num.diff<-dim(alignment[c(taxa1,taxa2),ind])[2]		# if there's 0 or 1 diff, this is going to give error
   for(i in 1:num.diff){
      aa.1<-unlist(strsplit(alignment[taxa1,ind[i]],""))
      aa.2<-unlist(strsplit(alignment[taxa2,ind[i]],""))
      diff<-which(aa.1!=aa.2)
      if(length(diff)==1){
         if(rbinom(1,1,perc)==1){
            free<-nucl[nucl!=aa.1[diff] & nucl!=aa.2[diff]]	# these are the ones that aren't used by the AAs
            test.split1<-aa.1					# or aa.2, doesn't matter
            test.split1[diff]<-free[1]		
            test.split2<-aa.1					# or aa.2, doesn't matter
            test.split2[diff]<-free[2]		

            test.comb1<-paste(test.split1[1],test.split1[2],test.split1[3],sep="")
            test.comb2<-paste(test.split2[1],test.split2[2],test.split2[3],sep="")

            both<-c(test.comb1,test.comb2)
            nonstop<-both[which(!is.element(both,stop.codons))]
            nonstop<-nonstop[which(!is.element(nonstop,alignment[-c(taxa1,taxa2),ind[i]]))]

            if(length(nonstop)>=1){
               replace<-sample(nonstop,1)	# if both are not stop codons, pick one

               align.taxa.conv[taxa1,ind[i]]<-replace
               align.taxa.conv[taxa2,ind[i]]<-replace

               others<-align.taxa.conv[-c(taxa1,taxa2),ind[i]]
               distinct<-length(unique(others))

               if(distinct>1){		# if they're all the same, don't do anything
                  if(distinct>2){	# if they're all different, just pick one at random
                     align.taxa.conv[-c(taxa1,taxa2),ind[i]]<-sample(others,1)	# i think... ???
                  } else {		# if there's a majority, make them all that
                     major<-which(colSums(sapply(others,"==",align.taxa.conv[,ind[i]]))==2)
                     align.taxa.conv[-c(taxa1,taxa2),ind[i]]<-others[major][1]
                  }
               }
            }
         }
      }
   }
   return(align.taxa.conv)
}
