make.mj.align <-
function(alignment){
   n.taxa<-nrow(alignment)
   n.nuc<-ncol(alignment)

   not.in.alphabet <- ((alignment != "A")*(alignment != "G")*(alignment != "C")*(alignment != "T")*
                      (alignment != "a")*(alignment != "g")*(alignment != "c")*(alignment != "t"))
   not.in.alphabet <- not.in.alphabet == 1

   align.table = matrix(0, nrow = n.taxa, ncol = n.nuc)
   rownames(align.table) = rownames(alignment)
    ##print(not.in.alphabet)

   align.table[not.in.alphabet] = 0
   align.table[alignment == "A" | alignment == "a"] = 1
   align.table[alignment == "G" | alignment == "g"] = 2
   align.table[alignment == "C" | alignment == "c"] = 3
   align.table[alignment == "T" | alignment == "t"] = 4
   return(align.table)
}


