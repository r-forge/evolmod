\name{plottree.db}
\alias{plottree.db}
\title{Plots the most likely phylogenetic trees}
\description{
  This function plots the most likely phylogenetic trees from a DualBrothers object.
}
\usage{plottree.db(x,numplot=4,makepic="FALSE",ext="png",type="unrooted",colors=NULL,threshold=0)}
\arguments{
  \item{x}{a DualBrothers object}
  \item{numplot}{a positive integer, determines the maximum number of trees to plot}
  \item{makepic}{a string, either "FALSE" or "TRUE", if "TRUE" the plot is outputted into a .png file in the current working directory}
  \item{ext}{a string, either "png" or "pdf", when makepic is TRUE determines what type of file to produce}
  \item{type}{a string specifying the type of phylogeny to be drawn; it must be one of "phylogram", "cladogram", "fan", "unrooted", "radial" or any unambiguous abbreviation of these}
   \item{colors}{a vector of color names to determine the palette for plotting}
  \item{threshold}{posterior probability threshold, a number between 0 and 1, only those trees that attained a posterior probability greater than the number provided will be plotted}

}
\value{Tree plots}
