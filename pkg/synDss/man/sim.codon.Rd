\name{sim.codon}
\alias{sim.codon}
\title{Simulate a codon sequence alignment}
\description{
   This function uses PAML to simulate sequences under codon models.
}
\usage{
sim.codon(n.taxa, n.codons, reps, tree, template.dat, datfile, freqs = c(0.6, 0.3, 0.1), omega = c(0.1, 0.8, 3.2), kap = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.taxa}{Number of taxa}
  \item{n.codons}{Number of codons}
  \item{reps}{Number of replicate datasets to simulate}
  \item{tree}{Phylogenetic tree}
  \item{template.dat}{Name of the file to use as the template for PAML}
  \item{datfile}{}
  \item{freqs}{The p parameter}
  \item{omega}{The omega paramater}
  \item{kap}{The kappa parameter}
}

\value{
   A list of sequence alignments
}
\references{
Z. Yang. PAML 4: phylogenetic analysis by maximum likelihood. \emph{Molecular Biology and Evolution}, \bold{24}:1586-1591, 2007.
}
\author{Peter Chi \email{pchi@calpoly.edu}
}
