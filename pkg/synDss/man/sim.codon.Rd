\name{sim.codon}
\alias{sim.codon}
\title{Simulate a codon sequence alignment}
\description{
   This function uses PAML to simulate sequences under codon models.
}
\usage{
sim.codon(n.taxa, n.codons, reps, tree, freqs = c(0.6, 0.3, 0.1), omega = c(0.1, 0.8, 3.2), kap = 2, call="evolverNSsites 6")
}
\arguments{
  \item{n.taxa}{Number of taxa}
  \item{n.codons}{Number of codons}
  \item{reps}{Number of replicate datasets to simulate}
  \item{tree}{Phylogenetic tree, in Newick format.}
  \item{freqs}{The p parameter}
  \item{omega}{The omega paramater}
  \item{kap}{The kappa parameter}
  \item{call}{The call to the evolver executable, with option "6" by default to indicate a simulation of codon sequences}
}

\details{Note: the phylogenetic tree must be represented by a character string in the Newick format. If your
tree is an object of class \code{"phylo"}, you must first use the \code{write.tree} function from \pkg{ape} to convert
it.}

\value{
   A list of simulated codon alignments
}
\references{
Z. Yang. PAML 4: phylogenetic analysis by maximum likelihood. \emph{Molecular Biology and Evolution}, \bold{24}:1586-1591, 2007.
}
\author{Peter Chi \email{pchi@calpoly.edu}
}
