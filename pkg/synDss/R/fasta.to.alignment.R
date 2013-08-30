fasta.to.alignment <-
function(fasta.file.name){
   return(ape::as.alignment(read.fasta(data,forceDNAtolower=F)))
}

# NOT SURE IF DONE YET, ALSO NEED .Rd file

