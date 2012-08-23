\name{plot.db}
\alias{plot.db}
\title{Plot a dualbrothers object}
\description{
  This function creates two plots from a DualBrothers object. The first plots the posterior probability of each tree as a function of base pair position. The second plots the posterior probability of a break point as a function of base pair position.
}
\usage{plot.db(x,makepic="FALSE",ext="png",seetrees="FALSE",numplot=4)}
\arguments{
  \item{x}{a DualBrothers object}
  \item{makepic}{a string, either "FALSE" or "TRUE", if "TRUE" the two plots are outputted into a .png file in the current working directory}
  \item{ext}{a string, either "png" or "pdf", when makepic is TRUE determines what type of file to produce}
  \item{seetrees}{a string, either "FALSE" or "TRUE", if "TRUE" the second plot show the trees instead of the posterior probability of a breakpoint}
  \item{numplot}{an integer, determines the number of trees to display in the second plot (only used if seestrees is "TRUE")}
  \item{type}{a string specifying the type of phylogeny to be drawn; it must be one of "phylogram" (the default),
          "cladogram", "fan", "unrooted", "radial" or any unambiguous abbreviation of these (only used if seestrees is "TRUE")}
}
\value{Two plots}
