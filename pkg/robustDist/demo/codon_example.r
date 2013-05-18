library(robustDist)

## set the seed
set.seed(39487348)

## read-in HIV codon-based multiple sequence alignment

my.align = read.phylip(file=system.file("sequences/pol.phy", package =
"robustDist"), format = "sequential", as.character = TRUE)

## convert nucleotides to codons

my.codon.align = dna.to.codon.align(my.align)

## define synonymous and nonsynonymous register matrices

syn.matrix = regist.synonym()
nonsyn.matrix = regist.nonsynonym()

## set the number of bootstrap simulations, this of course should be higher in a real analysis

sim.num = 10

## prepair two lists to store synonymous and nonsynonymous trees

syn.trees = list(sim.num)
nonsyn.trees = list(sim.num)

## run bootstrap simulations

for (i in 1:sim.num){

  ## resample codons in both nucleotide and codon alignments
  ## WARNING: we are resampling codon columns, not nucleotide columns
  
  resample.pair = resample.codon(my.align, my.codon.align)

  ## estimate synonymous and nonsynonymous distances

  syn.dist = robust.codon.dist.f84(resample.pair[[1]], syn.matrix, partition=1,resample.pair[[2]])
  nonsyn.dist = robust.codon.dist.f84(resample.pair[[1]], nonsyn.matrix, partition=1,resample.pair[[2]])
  

  syn.trees[[i]] = nj(syn.dist)
  nonsyn.trees[[i]] = nj(nonsyn.dist)
  
}

## uncommenting the next two lines will allow you to write the bootstrap trees into
## two nexus files and then study them using your favorite tree manipulation software
## we however are going to proceed with ape capabilities of tree manipulation

#write.nexus(syn.trees, file = "syn.nex", translate = TRUE, original.data = FALSE)
#write.nexus(nonsyn.trees, file = "nonsyn.nex", translate = TRUE, original.data = FALSE)

## build majority consesus tree

syn.cons = consensus(syn.trees,p=0.5)
nonsyn.cons = consensus(nonsyn.trees,p=0.5)

## compute bootstrap support of the consesus trees bipartitions and attach these
## probabilities to the nodes of the consesus trees

syn.cons$node.label = round(100*prop.clades(syn.cons, syn.trees)/length(syn.trees))
nonsyn.cons$node.label = round(100*prop.clades(nonsyn.cons, nonsyn.trees)/length(nonsyn.trees))

boot.color = "lightblue"

par(mfrow=c(1,2), mar=c(5,0,4,0)+0.1)

plot(syn.cons, type = "u", show.node.label=FALSE,edge.width=3,cex=0.8,
     main = "Synonymous Distance Tree", cex.main = 1.2)

nodelabels(text=syn.cons$node.label, cex = 0.8, frame = "c", bg = boot.color, font = 2)


plot(nonsyn.cons, type = "u", show.node.label=FALSE,edge.width=3,cex = 0.8,
     main = "Nonsynonymous Distance Tree", cex.main = 1.2)

nodelabels(text=nonsyn.cons$node.label, cex = 0.8, frame = "c", bg = boot.color,
           font = 2)

