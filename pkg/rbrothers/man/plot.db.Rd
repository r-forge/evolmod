\name{plot.db}
\alias{plot.db}
\title{Plot a dualbrothers object}
\description{
  This function creates two plots from a DualBrothers object. The first plots the posterior probability of each tree as a function of base pair position. The second plots the posterior probability of a break point as a function of base pair position.
}
\usage{\method{plot}{db}(x,seetrees="FALSE",numplot=4,type="unrooted",colors=NULL, ...)}
\arguments{
  \item{x}{a DualBrothers object}
  \item{makepic}{a string, either "FALSE" or "TRUE", if "TRUE" the two plots are outputted into a .png file in the current working directory}
  \item{ext}{a string, either "png" or "pdf", when makepic is TRUE determines what type of file to produce}
  \item{seetrees}{a string, either "FALSE" or "TRUE", if "TRUE" the second plot show the trees instead of the posterior probability of a breakpoint}
  \item{numplot}{an integer, determines the number of trees to display in the second plot (only used if seestrees is "TRUE")}
  \item{type}{a string specifying the type of phylogeny to be drawn; it must be one of "phylogram", "cladogram", "fan", "unrooted", "radial" or any unambiguous abbreviation of these (only used if seestrees is "TRUE")}
  \item{linewidth}{a number, determines the thickness of the lines on the posterior probability plot}
  \item{colors}{a vector of color names to determine the palette for plotting}
  \item{threshold}{posterior probability threshold, a number between 0 and 1, only those trees that attained a posterior probability greater than the number provided will be plotted}
  \item{...}{optional arguments}
}
\value{Two posterior probability plots}
