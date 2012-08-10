\name{plot.db}
\alias{plot.db}
\title{Plot a dualbrothers object}
\description{
  This function creates two plots from a DualBrothers object. The first plots the posterior probability of each tree as a function of base pair position. The second plots the posterior probability of a break point as a function of base pair position.
}
\usage{plot.db(x,makepng="FALSE")}
\arguments{
  \item{x}{a DualBrothers object}
  \item{makepng}{a string, either "FALSE" or "TRUE", if "TRUE" the two plots are outputted into a .png file in the current working directory}
}
\value{Two plots}
