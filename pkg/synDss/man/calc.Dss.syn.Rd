\name{calc.Dss.syn}
\alias{calc.Dss.syn}
\title{Calculate the Dss statistic using synonymous distances}
\description{
   This function calculates the Dss statistic using sliding windows and synonymous distance
}
\usage{
calc.Dss.syn(data, l, m, syn.matrix, exhaustive=FALSE)
}

\arguments{
  \item{data}{The sequence alignment, as a matrix of A,C,T,G characters}
  \item{l}{Step size for the sliding windows}
  \item{m}{Half-window size (total window size will be 2m)}
  \item{syn.matrix}{An indicator matrix which specifies the transitions of interest. For synonymous distances,
each entry will be a 1 for synonymous changes, and 0 for non-synonymous changes. See \code{\link{regist.synonym}}.}
  \item{exhaustive}{logical; if TRUE, an exhaustive search over all possible topologies is performed at each step.}
}
\details{
  Unless \code{exhaustive=TRUE}, this function uses \code{\link{bionj}} to quickly obtain the neighbor-joining tree before finding
optimal branch lengths using least squares for each half-window; the branch length estimation is restricted to the topology obtained by the
neighbor-joining algorithm.
}
\value{
  A vector of numeric values, representing the Dss statistic at each window location.
}
\references{
G. McGuire, F. Wright, and M.J. Prentice. A graphical method for detecting recombination
in phylogenetic data sets. \emph{Molecular Biology and Evolution}, \bold{14}:1125-1131, 1997.

N. Saitou and M. Nei. The neighbor-joining method: A new method for reconstructing
phylogenetic trees. \emph{Molecular Biology and Evolution}, \bold{4}:406-425, 1987.

J.D. O'Brien, V.N. Minin, and M.A. Suchard. Learning to count: robust estimates for
labeled distances between molecular sequences. \emph{Molecular Biology and Evolution}, \bold{26}:801-814, 2009.
}
\author{Peter Chi \email{pchi@calpoly.edu}
}
\seealso{
  \code{\link{calc.Dss}}
}
\keyword{recombination}
\keyword{Dss}
\keyword{phylogenetics}
\keyword{synonymous}
