optim.phylo.ls.all <-
function (D, set.neg.to.zero = TRUE, fixed = FALSE, 
    tol = 1e-10, collapse = TRUE) 
{
    if (class(D) == "dist") 
        D <- as.matrix(D)
    n <- nrow(D)

    all.trees <- allTrees(n,tip.label=row.names(as.matrix(D)))
    allQ <- vector()
    bestQ <- Inf
    for (i in 1:length(all.trees)) {
        all.trees[[i]]$edge.length<-rep(dim(all.trees[[i]]$edge)[1],1)
        all.trees[[i]] <- nnls.tree(D,all.trees[[i]],trace=0)	# this used to be ls.tree
        allQ[i] <- attr(all.trees[[i]], "RSS")
    }
    best<-which(allQ==min(allQ))
    if(length(best)>1){			# temporary fix?
        best<-best[1]			# if trees tie, just pick the first one
    }
    return(all.trees[[best]])
}
