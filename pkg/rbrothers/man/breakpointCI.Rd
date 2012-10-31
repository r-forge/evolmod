\name{breakpointCI}
\alias{breakpointCI}
\title{95\% breakpoint credible interval}
\description{
Returns the 95\% credible interval for a single break point in between the user provided lower and upper bounds (nucleotide position) for a DualBrothers object, x.
}
\usage{breakpointCI(x,i,j)}
\arguments{
  \item{x}{a dualbrothers object}
  \item{i}{an integer indicating the nucleotide position of the left edge of an interval containing one breakpoint}
  \item{j}{an integer indicating the nucleotide position of the right edge of an interval containing one breakpoint}
}
\value{The 95\% credible interval for a single breakpoint in the interval provided.}
