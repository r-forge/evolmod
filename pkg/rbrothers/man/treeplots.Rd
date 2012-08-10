\name{treeplots}
\alias{treeplots}
\title{DualBrothers Tree Plot Maker}
\description{
  This function creates a plot of the top trees found by DualBrothers.
}
\usage{treeplots("basename",numplot=12,visualize="FALSE")}
\arguments{
  \item{basename}{a string in quotes ("basename") that requires the basename.topoprob, basename.profile, and basename.tree files to be in the current working directory}
 \item{numplot}{an integer that determines the maximum number of trees to plot, defaults to 12}
 \item{visualize}{either TRUE or FALSE (default is FALSE); if TRUE, the tree plot will be shown}
}
\value{One .png plot is made in the current working directory.}
