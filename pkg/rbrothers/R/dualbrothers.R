"dualbrothers"<- function(seed,alignment,format="interleaved",window.size=NULL,basename=alignment,length=2100000,burnin=100000,subsample=200,par_lambda=5,top_lambda=0.693,window_length=10,C=0.3,sigma_alpha=0.75,sigma_mu=0.75,subst_hyper_mean=-1,subst_hyper_variance=-1,diver_hyper_mean=-1,diver_hyper_variance=-1,top_breaks=1,par_breaks=1,step.size=10,boot=0){
  
   options(scipen = 100)

   write(paste("length:",length), file = paste(basename,".cmdfile",sep=""))
   cat(paste("burnin:",burnin,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   cat(paste("subsample:",subsample,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   cat(paste("par_lambda:",par_lambda,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   cat(paste("top_lambda:",top_lambda,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   cat(paste("window_length:",window_length,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   cat(paste("C:",C,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   cat(paste("sigma_alpha:",sigma_alpha,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   cat(paste("sigma_mu:",sigma_mu,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE) 
   cat(paste("top_breaks:",top_breaks,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   cat(paste("par_breaks:",par_breaks,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   if(subst_hyper_mean!=-1&subst_hyper_variance!=-1&diver_hyper_mean!=-1&diver_hyper_variance!=-1){
    cat(paste("subst_hyper_mean:",subst_hyper_mean,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
    cat(paste("subst_hyper_variance:",subst_hyper_variance,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE) 
    cat(paste("diver_hyper_mean:",diver_hyper_mean,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
    cat(paste("diver_hyper_variance:",diver_hyper_variance,"\n"), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   }
   if(subst_hyper_mean!=-1|subst_hyper_variance!=-1|diver_hyper_mean!=-1|diver_hyper_variance!=-1){
    if(subst_hyper_mean==-1|subst_hyper_variance==-1|diver_hyper_mean==-1|diver_hyper_variance==-1){
     print("All four hyper variables must be provided. At least one is missing so they will all be estimated.")
    }
   }

library(ape)

write.dna<-function (x, file, format = "interleaved", append = FALSE, nbcol = 6,colsep = " ", colw = 10, indent = NULL, blocksep = 1) 
{
    format <- match.arg(format, c("interleaved", "sequential", 
        "fasta"))
    phylip <- if (format %in% c("interleaved", "sequential")) 
        TRUE
    else FALSE
    if (inherits(x, "DNAbin")) 
        x <- as.character(x)
    aligned <- TRUE
    if (is.matrix(x)) {
        N <- dim(x)
        S <- N[2]
        N <- N[1]
        xx <- vector("list", N)
        for (i in 1:N) xx[[i]] <- x[i, ]
        names(xx) <- rownames(x)
        x <- xx
        rm(xx)
    }
    else {
        N <- length(x)
        S <- unique(unlist(lapply(x, length)))
        if (length(S) > 1) 
            aligned <- FALSE
    }
    if (is.null(names(x))) 
        names(x) <- as.character(1:N)
    if (is.null(indent)) 
        indent <- if (phylip) 
            10
        else 0
    if (is.numeric(indent)) 
        indent <- paste(rep(" ", indent), collapse = "")
    if (format == "interleaved") {
        blocksep <- paste(rep("\n", blocksep), collapse = "")
        if (nbcol < 0) 
            format <- "sequential"
    }
    zz <- if (append) 
        file(file, "a")
    else file(file, "w")
    on.exit(close(zz))
    if (phylip) {
        if (!aligned) 
            stop("sequences must have the same length for\n interleaved or sequential format.")
        cat(N, " ", S, "\n", sep = "", file = zz)
        if (nbcol < 0) {
            nb.block <- 1
            nbcol <- totalcol <- ceiling(S/colw)
        }
        else {
            nb.block <- ceiling(S/(colw * nbcol))
            totalcol <- ceiling(S/colw)
        }
        SEQ <- matrix("", N, totalcol)
        for (i in 1:N) {
            X <- paste(toupper(x[[i]]), collapse = "")
            for (j in 1:totalcol) SEQ[i, j] <- substr(X, 1 + 
                (j - 1) * colw, colw + (j - 1) * colw)
        }
        max.nc <- max(nchar(names(x)))
        if (max.nc < 11) 
            max.nc <- 9
        fmt <- paste("%-", max.nc + 1, "s", sep = "")
        names(x) <- sprintf(fmt, names(x))
    }
    if (format == "interleaved") {
        colsel <- if (nb.block == 1) 
            1:totalcol
        else 1:nbcol
        for (i in 1:N) {
            cat(paste(names(x)[i]," ",sep=""), file = zz)
            cat(SEQ[i, colsel], sep = colsep, file = zz)
            cat("\n", file = zz)
        }
        if (nb.block > 1) {
            for (k in 2:nb.block) {
                cat(blocksep, file = zz)
                endcolsel <- if (k == nb.block) 
                  totalcol
                else nbcol + (k - 1) * nbcol
                for (i in 1:N) {
                  cat(indent, file = zz)
                  cat(SEQ[i, (1 + (k - 1) * nbcol):endcolsel], 
                    sep = colsep, file = zz)
                  cat("\n", file = zz)
                }
            }
        }
    }
    if (format == "sequential") {
        if (nb.block == 1) {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, ], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        }
        else {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
                cat("\n", file = zz)
                for (k in 2:nb.block) {
                  endcolsel <- if (k == nb.block) 
                    totalcol
                  else nbcol + (k - 1) * nbcol
                  cat(indent, file = zz)
                  cat(SEQ[i, (1 + (k - 1) * nbcol):endcolsel], 
                    sep = colsep, file = zz)
                  cat("\n", file = zz)
                }
            }
        }
    }
    if (format == "fasta") {
        for (i in 1:N) {
            cat(">", names(x)[i], file = zz)
            cat("\n", file = zz)
            X <- paste(x[[i]], collapse = "")
            S <- length(x[[i]])
            totalcol <- ceiling(S/colw)
            if (nbcol < 0) 
                nbcol <- totalcol
            nb.lines <- ceiling(totalcol/nbcol)
            SEQ <- character(totalcol)
            for (j in 1:totalcol) SEQ[j] <- substr(X, 1 + (j - 
                1) * colw, colw + (j - 1) * colw)
            for (k in 1:nb.lines) {
                endsel <- if (k == nb.lines) 
                  length(SEQ)
                else nbcol + (k - 1) * nbcol
                cat(indent, file = zz)
                cat(SEQ[(1 + (k - 1) * nbcol):endsel], sep = colsep, 
                  file = zz)
                cat("\n", file = zz)
            }
        }
    }
}




if(format=="fasta"){
small.align = read.dna(paste(alignment,".fasta",sep=""), format = "fasta")
write.dna(small.align,file=paste(basename,".phy",sep=""))
}
if(format=="interleaved"){
small.align = read.dna(paste(alignment,".phy",sep=""), format = "interleaved")
write.dna(small.align,file=paste(basename,".phy",sep=""))
}

align.length = dim(small.align)[2]
taxa.num = dim(small.align)[1]
taxa.names = rownames(small.align)
rownames(small.align) = c(1:taxa.num)-1
sw.trees = list(1)
my.counter = 1

if(taxa.num>6){

 if(window.size==NULL){
  print("'window.size' is missing")
  return()
 }

for(ws in 1:length(window.size)){
 window.start = 1
 window.end = window.start + window.size[ws] -1
 while (window.end < align.length){
  temp.tree = bionj(dist.dna(small.align[,window.start:window.end]))
  temp.tree$edge.length = NULL
  sw.trees[[my.counter]] = temp.tree
  window.start = window.start + step.size
  window.end = window.end + step.size
  my.counter = my.counter + 1
 }
}

#### bootstrap the dna
if(boot>0){
 for(k in 1:boot){
  for(ws in 1:length(window.size)){
   window.start = 1
   window.end = window.start + window.size[ws] -1
   while (window.end < align.length){
    temp.tree = bionj(dist.dna(small.align[,sample(window.start:window.end,replace=TRUE)]))
    temp.tree$edge.length = NULL
    sw.trees[[my.counter]] = temp.tree
    window.start = window.start + step.size
    window.end = window.end + step.size
    my.counter = my.counter + 1
   }
  }
 }
}

cleaned.sw.trees = unique(sw.trees)
}

if(taxa.num==6){cleaned.sw.trees = read.tree(paste(system.file(package="rbrothers"),"/extdata/smalltrees/six-input.tre",sep=""))}
if(taxa.num==5){cleaned.sw.trees = read.tree(paste(system.file(package="rbrothers"),"/extdata/smalltrees/five-input.tre",sep=""))}
if(taxa.num==4){cleaned.sw.trees = read.tree(paste(system.file(package="rbrothers"),"/extdata/smalltrees/four-input.tre",sep=""))}
if(taxa.num==3){cleaned.sw.trees = read.tree(paste(system.file(package="rbrothers"),"/extdata/smalltrees/three-input.tre",sep=""))}

tree.num = length(cleaned.sw.trees)

for (i in 1:tree.num){
  my.tree = root(cleaned.sw.trees[[i]], outgroup="0", r=TRUE)
  write.tree(my.tree, file=paste(basename,"-input.tre",sep=""), append=TRUE)
}

temp.tree = read.tree(paste(basename,"-input.tre",sep=""))
unique.trees = unique(temp.tree)
write.tree(unique.trees, file=paste(basename,"-input.tre",sep=""))


   cat(paste("tree_file: ",basename,"-input.tre",sep=""), file = paste(basename,".cmdfile",sep=""),append=TRUE)
   dbs <- .jnew("DualBrothersForR")
   out <- .jcall(dbs, "V", "main", c(seed,paste(basename,".cmdfile",sep=""),paste(basename,".phy",sep=""),paste(basename,".post",sep=""))) 
   #system(paste("java DualBrothers ",seed," ",basename,".cmdfile ",basename, ".phy ",basename,".post",sep=""))
   topologyprofile(basename)
   epprofile(basename)

   
probs<-read.table(paste(basename,".topoprob",sep=""))
profile<-read.table(paste(basename,".profile",sep=""))

#### creates a vector of all the tip names culled from the .phy file
P<-as.integer(scan(paste(basename,".phy",sep=""),what='raw',nlines=1)[1])
names<-NULL
for(i in 1:P){
names[i]<-(scan(paste(basename,".phy",sep=""),what='raw',nlines=1,skip=i)[1])
}

breaks<-findbreakpoints(basename)

#### read in the trees
treenum<-dim(probs)[2]-1
trees<-read.tree(paste(basename,".tree",sep=""),tree.names=seq(1:treenum),keep.multi=TRUE)

#### relabel the tips with the names from the .phy file
for(i in 1:treenum){
 for(j in 1:P){
  trees[[i]]$tip.label[j]<-names[as.integer(trees[[i]]$tip.label[j])+1]
 }
}

####
#### read the .post file to find the number of topology breakpoints and log likelihoods
####

lineiloglik<-function(linei){return(sum(as.numeric(linei[5+8*(0:(as.numeric(linei[2])-1))])))}

lineinumbrkpnts<-function(linei){
 numchgpnts<-as.numeric(linei[2])
 ts1<-linei[4]
 returnme<-0
 for(i in 2:numchgpnts){
  if(!linei[4+8*(i-2)]==linei[4+8*(i-1)]){returnme<-returnme+1}
 }
 return(returnme)
}

loglikes<-NULL
numbrkpnts<-NULL
ll <- readLines(paste(basename,".post",sep=""))

for(i in 1:length(ll)){
 ii = unlist(strsplit(ll[[i]]," "))
 ii = ii[which(sapply(ii,nchar,USE.NAMES=FALSE)!=0)]
 loglikes[i]<-lineiloglik(ii)
 numbrkpnts[i]<-lineinumbrkpnts(ii)
}

####
####
####

tbr<-list(basename,probs,profile,P,breaks,trees,par_lambda,top_lambda,sigma_alpha,sigma_mu,top_breaks,par_breaks,loglikes,numbrkpnts)
class(tbr)<-"db"
names(tbr)<-c("basename","TopologyProfile","EPProfile","numberofsequences","breakpoints","trees","par_lambda","top_lambda","sigma_alpha","sigma_mu","top_breaks","par_breaks","MCMC_log_likelihood","MCMC_number_of_break_points")
return(tbr)

}






