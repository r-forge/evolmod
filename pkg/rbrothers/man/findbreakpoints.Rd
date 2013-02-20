\name{findbreakpoints}
\alias{findbreakpoints}
\title{FindBreakPoints}
\description{
   This function is used to find the sites where the tree with the largest posterior probability changes. A basename.topoprob file is needed. This file is created by the TopologyProfile function. 
}
\usage{findbreakpoints("basename")}
\arguments{
  \item{basename}{a string in quotes that requires the basename.topoprob file to be in the current working directory}
}
\details{Findbreakpoints returns a vector containing the site numbers where the posterior probability of the largest tree changes.) 
}


