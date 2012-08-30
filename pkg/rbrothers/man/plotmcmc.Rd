\name{plotmcmc}
\alias{plotmcmc}
\title{MCMC trace plot of log likelihoods}
\description{
  This function uses the coda package to produce a trace plot of the log likelihoods in the mcmc run. It can also produce an autocorrelation plot of the log likelihoods.
}
\usage{plotmcmc(x,mode="trace")}
\arguments{
  \item{x}{a dualbrothers object}
  \item{mode}{a string, either "trace" or "auto" determining which type of plot to create}
}
\value{A trace plot or an autocorrelation plot of the MCMC log likelihoods.}
