\name{calc.null.pboot.syn}
\alias{calc.null.pboot.syn}
\title{Calculate the null distribution of the Synonymous Dss statistic using parametric bootstrap}
\description{
   This function estimates a null distribution of the Dss statistic using the overall
least squares tree as the evolutionary history of the sequences under the null hypothesis. The null
distribution is calculated using the Synonyous Dss statistic across the alignment.
}
\usage{
calc.null.pboot.syn(data, B, l, m, syn.matrix, exhaustive=F, codeml="codeml", evolver="evolverNSsites 6")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The sequence alignment, as a matrix of A,C,T,G characters}
  \item{B}{Number of bootstrap replicates to use}
  \item{l}{Step size for the sliding windows}
  \item{m}{Half-window size (total window size will be 2m)}
  \item{syn.matrix}{An indicator matrix which specifies the transitions of interest. For synonymous distances,
each entry will be a 1 for synonymous changes, and 0 for non-synonymous changes. See \code{\link{regist.synonym}}.} 
  \item{exhaustive}{Logical: If TRUE, then an exhaustive search over all possible topologies is performed. Not allowed
for datasets of greater than 6 sequences due to memory/time constraints.}
  \item{codeml}{Call to codeml executable; can be changed if needed}
  \item{evolver}{Call to evolver executable; can be changed if needed}
}
\value{
  A vector of numeric values, representing the maximum Dss statistic at each bootstrap replicate.
}
\references{
G. McGuire, F. Wright, and M.J. Prentice. A graphical method for detecting recombination
in phylogenetic data sets. \emph{Molecular Biology and Evolution}, \bold{14}:1125-1131, 1997.
}
\author{Peter Chi \email{pchi@calpoly.edu}
}

\seealso{
\code{\link{calc.Dss.syn}}, \code{\link{calc.null.pboot}}
}

\keyword{recombination}
\keyword{Dss}
\keyword{phylogenetics}
