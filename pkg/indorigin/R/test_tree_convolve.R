treeConvolveTest= function(my.tree, my.data, rate.01, rate.10,
		root.dist.prob.0, n.max) {
	if (!("phylo" %in% class(my.tree)))
		stop("Error: object \"my.tree\" is not of class \"phylo\"")
	
	if (is.null(my.tree$edge.length))
		stop("Error: tree \" my.tree\" must have branch lengths.")
	
	if (!is.rooted(my.tree))
	stop("Error: The input tree must be rooted")
  
	## reorder data on tips to match the order of the my.tree phylo object
	my.tree = reorder(my.tree, order = "postorder")
	if (!is.null(names(my.data))) {
		if (!any(is.na(match(names(my.data), my.tree$tip.label)))) {
			my.data = my.data[my.tree$tip.label]
		} else {
			warning('the names of argument "my.data" and the names of the tip labels
							did not match: the former were ignored in the analysis.')
		}
	}

	res = .Call("treeConvolve", my.tree$edge, my.tree$edge.length,
			my.data, rate.01, rate.10, n.max, root.dist.prob.0,
			PACKAGE = "indorigin")
	colnames(res) = c("prior", "posterior")
	return(res)
}