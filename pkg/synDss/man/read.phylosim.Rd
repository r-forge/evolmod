\name{read.phylosim}
\alias{read.phylosim}
\title{Convert alignment to be readable for synonymous distance calculations}
\description{
  The synonymous distance function in the robustDist package uses a 1,2,3,4 coding
for nucleotides instead of A,C,G,T. This function makes the conversion.}
\usage{
read.phylosim(alignment)
}
\arguments{
  \item{alignment}{A sequence alignment in the format of a matrix of A,C,G,T}
}
\value{
  Matrix of 1,2,3,4
}

