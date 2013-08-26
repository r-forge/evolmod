\name{codeml.M3}
\alias{codeml.M3}
\title{Estimate parameters from the M3 codon substitution model}
\description{
  This function takes a sequence alignment dataset and an assumed tree,
and infers the M3 model parameters kappa, p and omega (p and omega are vectors). Internal
function, used within parametric bootstrap.}
\usage{
codeml.M3(id, data, ctl.name, tree)
}
\arguments{
  \item{id}{Given from parametric bootstrap function}
  \item{data}{Sequence alignment dataset}
  \item{ctl.name}{Name of the control file for paml}
  \item{tree}{Phylogenetic tree, the assumed evolutionary history.}
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
