\name{plottree.db}
\alias{plottree.db}
\title{Plots the most likely phylogenetic trees}
\description{
  This function plots the most likely phylogenetic trees from a DualBrothers object.
}
\usage{plottree.db(x,numplot=12,makepic="FALSE",ext="png")}
\arguments{
  \item{x}{a DualBrothers object}
  \item{numplot}{a positive integer, determines the maximum number of trees to plot}
  \item{makepic}{a string, either "FALSE" or "TRUE", if "TRUE" the plot is outputted into a .png file in the current working directory}
  \item{ext}{a string, either "png" or "pdf", when makepic is TRUE determines what type of file to produce}
}
\value{Tree plots}
