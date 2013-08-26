\name{optim.phylo.ls.all}
\alias{optim.phylo.ls.all}
\title{Least squares phylogenetic estimation with an exhaustive search over all possible topologies}
\description{
   This function obtains the phylogeny which has the minimum value of the least squares
criterion, among all possible phylogenies. An exhaustive search over all possible
topologies is performed, to ensure that a global minimum is reached in the optimization routine.
}
\usage{
optim.phylo.ls.all(D, set.neg.to.zero = TRUE, fixed = FALSE, tol = 1e-10, collapse = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{Distance matrix}
  \item{set.neg.to.zero}{logical; if TRUE, set negative values equal to 0.}
  \item{fixed}{logical; if TRUE, use a fixed topology for branch length estimation}
  \item{tol}{Tolerance}
  \item{collapse}{logical; if TRUE, then collapse branches with zero length}
}
\details{
  This function is a modification of \code{optim.phylo.ls} from \pkg{phytools}
}
\value{
   An object of class \code{"phylo"} that is the least-squares tree with branch lengths; 
also returns the sum of squares in \code{attr(tree,"Q-score")}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. \emph{3}, 217-223. doi:10.1111/j.2041-210X.2011.00169.x
}
\author{
Peter Chi
}

