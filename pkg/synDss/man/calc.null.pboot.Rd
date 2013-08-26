\name{calc.null.pboot}
\alias{calc.null.pboot}
\title{Calculate the null distribution of the Dss statistic using parametric bootstrap}
\description{
   This function takes the data, and estimates a null distribution of the Dss statistic
using the overall least squares tree as the evolutionary history of the sequences under the null
hypothesis.
}
\usage{
calc.null.pboot(data, B, l, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The sequence alignment, as a matrix of A,C,T,G characters}
  \item{B}{Number of bootstrap replicates to use}
  \item{l}{Step size for the sliding windows (should be the same as what was used for calculation of the Dss statistic)}
  \item{m}{Half-window size (total window size will be 2m)}
}
\details{
  This function currently uses an exhaustive search over all possible topologies.
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
\code{\link{calc.Dss.all}}, \code{\link{calc.null.pboot.syn}}
}

\keyword{recombination}
\keyword{Dss}
\keyword{phylogenetics}
