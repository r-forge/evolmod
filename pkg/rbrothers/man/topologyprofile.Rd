\name{topologyprofile}
\alias{topologyprofile}
\title{TopologyProfile }
\description{
   This function is used to ease the interpretation of the outputfile, basename.post, created by dualbrothers. TopologyProfile creates the files basename.topoprob and basename.tree. 
}
\usage{topologyprofile("basename")}
\arguments{
  \item{basename}{A string in quotes that requires the basename.post file to be in the current working directory. The basename of an rbrothers object, x, can be extracted via x$basename.}
}
\details{TopologyProfile creates two files that are more readable than the output file created by DualBrothers, basename.post. Basename.topoprob encodes the posterior probabilities at each site in the sequence (the rows) for each tree (the columns). Each column corresponds to its respective tree in the basename.tree file. (The first column in the basename.topoprob file corresponds to the first tree in the basename.tree file.) 
}
\value{TopologyProfile creates a basename.topoprob file and a basename.tree file in the current working directory.}

