\name{plottree.db}
\alias{plottree.db}
\title{Plots the most likely phylogenetic trees}
\description{
  This function plots the most likely phylogenetic trees from a DualBrothers object.
}
\usage{plottree.db(x,numplot=4,makepic="FALSE",ext="png",type="unrooted")}
\arguments{
  \item{x}{a DualBrothers object}
  \item{numplot}{a positive integer, determines the maximum number of trees to plot}
  \item{makepic}{a string, either "FALSE" or "TRUE", if "TRUE" the plot is outputted into a .png file in the current working directory}
  \item{ext}{a string, either "png" or "pdf", when makepic is TRUE determines what type of file to produce}
  \item{type}{a string specifying the type of phylogeny to be drawn; it must be one of "phylogram", "cladogram", "fan", "unrooted", "radial" or any unambiguous abbreviation of these}
}
\value{Tree plots}
