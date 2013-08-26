\name{conv.evol}
\alias{conv.evol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate convergent evolution}
\description{
  This function simulates convergent evolution between two taxa from a DNA sequence alignment.
}
\usage{
conv.evol(alignment, perc, taxa1, taxa2, map)
}

\arguments{
  \item{alignment}{The DNA sequence alignment, as a matrix of A,C,T,G characters}
  \item{perc}{Percentage of possible one-step codon changes into identical amino acids to make.}
  \item{taxa1}{A number indicating the row of the first taxa to use}
  \item{taxa2}{A number indicating the row of the second taxa to use}
  \item{map}{Codon map}
}

