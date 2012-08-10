\name{epprofile}
\alias{epprofile}
\title{EPProfile}
\description{
  This function is used to ease the interpretation of the output file, basename.post, created by dualbrothers. EPProfile creates the file basename.profile.
}
\usage{epprofile("basename")}
\arguments{
  \item{basename}{a string in quotes ("basename") that requires the basename.post file and the basename.profile file to be in the current working directory}
}
\details{
The basename.profile file contains the posterior probability of the evolutionary parameters at each site. The format is as follows:
<site> <kappa-mean> <mu-mean> 
<kapp-q1.5%> <kappa-q50%> <kappa-q97.5%> <mu-q2.5%> 
<mu-q50%> <mu-q97.5%> <COP>.
}
\value{EPProfile creates a basename.profile file in the current working directory.}
