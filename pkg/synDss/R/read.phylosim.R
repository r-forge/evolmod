read.phylosim <-
function(alignment){
   n.taxa<-nrow(alignment)
   n.codons<-ncol(alignment)
   matrixed.align<-t(matrix(as.vector(unlist(apply(alignment,1,strsplit,""))),ncol=n.taxa))

   not.in.alphabet <- ((matrixed.align != "A")*(matrixed.align != "G")*(matrixed.align != "C")*(matrixed.align != "T")*
                      (matrixed.align != "a")*(matrixed.align != "g")*(matrixed.align != "c")*(matrixed.align != "t"))
   not.in.alphabet <- not.in.alphabet == 1

   align.table = matrix(0, nrow = n.taxa, ncol = n.codons*3)
   rownames(align.table) = rownames(alignment)
    ##print(not.in.alphabet)

   align.table[not.in.alphabet] = 0
   align.table[matrixed.align == "A" | matrixed.align == "a"] = 1
   align.table[matrixed.align == "G" | matrixed.align == "g"] = 2
   align.table[matrixed.align == "C" | matrixed.align == "c"] = 3
   align.table[matrixed.align == "T" | matrixed.align == "t"] = 4
   return(align.table)
}
