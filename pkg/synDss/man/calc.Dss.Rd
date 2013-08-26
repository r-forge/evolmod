\name{calc.Dss}
\alias{calc.Dss}
\title{Calculate the Dss statistic
}
\description{
   This function calculates the Dss statistic using sliding windows.
}
\usage{
calc.Dss(data, l, m)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{data}{The sequence alignment, as a matrix of A,C,T,G characters}
  \item{l}{Step size for the sliding windows}
  \item{m}{Half-window size (total window size will be $2 \times m$)}
}
\details{
  In contrast to calc.Dss.all, this function uses \code{\link{bionj}} to quickly obtain the neighbor-joining tree before finding
optimal branch lengths using least squares; the branch length optimization is restricted to the topology obtained by
neighbor-joining. This allows for a much faster implementation than calc.Dss.all
}
\value{
  A vector of numeric values, representing the Dss statistic at each window location.
}
\references{
G. McGuire, F. Wright, and M.J. Prentice. A graphical method for detecting recombination
in phylogenetic data sets. \emph{Molecular Biology and Evolution}, \bold{14}:1125-1131, 1997.

N. Saitou and M. Nei. The neighbor-joining method: A new method for reconstructing
phylogenetic trees. \emph{Molecular Biology and Evolution}, \bold{4}:406-425, 1987.

}
\author{Peter Chi \email{pchi@calpoly.edu}
}
\seealso{
  \code{\link{calc.Dss.all}}, \code{\link{calc.Dss.syn}}, \code{\link{calc.Dss.syn.all}}
}
\keyword{recombination}
\keyword{Dss}
\keyword{phylogenetics}

