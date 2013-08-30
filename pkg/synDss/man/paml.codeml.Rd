\name{paml.codeml}
\alias{paml.codeml}
\title{Estimate parameters from the M3 codon substitution model}
\description{
  This function takes a sequence alignment dataset and an assumed tree,
and infers the M3 model parameters kappa, p and omega (p and omega are vectors). Internal
function, used within parametric bootstrap.}
\usage{
paml.codeml(data, tree, call="codeml")
}
\arguments{
  \item{data}{Sequence alignment dataset, an object of class \code{"alignment"}}
  \item{tree}{Phylogenetic tree as an object of class \code{"phylo"}}
  \item{call}{Command to the codeml executable}
}

\value{
 Vector of kappa, p and omega
}
\references{
Z. Yang and R. Nielsen. Codon-substitution models for detecting molecular adaptation
at individual sites along specific lineages. \emph{Molecular Biology and Evolution},
\bold{19}:908-917, 2002.
}
\author{Peter Chi \email{pchi@calpoly.edu}
}
