\name{dualbrothers}
\alias{dualbrothers}
\title{Recombination Analysis via DualBrothers}
\description{
  DualBrothers generates a posterior sample based on the dual Multiple Change-Point model.
}
\usage{dualbrothers(seed,alignment,format="interleaved",window.size=NULL,basename=alignment,length=2100000,burnin=100000,subsample=200,par_lambda=5,top_lambda=0.693,window_length=10,C=0.3,sigma_alpha=0.75,sigma_mu=0.75,subst_hyper_mean=-1,subst_hyper_variance=-1,diver_hyper_mean=-1,diver_hyper_variance=-1,top_breaks=1,par_breaks=1, step.size = 10)}
\arguments{
  \item{seed}{integer used to seed the pseudo-random number generator (required)}
  \item{alignment}{a string in quotes ("alignment") that requires the alignment.fasta file or the alignment.phy file to be in the current working directory (required)}
  \item{format}{a string, either "fasta" or "interleaved"}
  \item{window.size}{a vector of positive integers indicating the widths of the windows (in base pairs) used to create the set of candidate trees (required if the number of taxa is greater than 6)}
  \item{basename}{a string in quotes ("basename") to be used for output files}
  \item{length}{a positive integer indicating the total number of samples to generate}
  \item{burnin}{a positive integer indicating the number of initial samples to discard}
  \item{subsample}{a positive number indicating the frequency at which samples are saved to the output file after the burnin}
  \item{par_lambda}{hyperprior parameter representing the prior mean number of substitution process change-points}
  \item{top_lambda}{hyperprior parameter representing the prior mean number of topology break-points}
  \item{window_length}{size of window around which existing change-points are randomly moved during the update of Peter Green's constant}
  \item{C}{Peter Green's constant. Adjust this constant to change the proportion of time that the sampler spends updating parameters when the number of dimensions is fixed.}
  \item{sigma_alpha}{spread of alphas when proposing new segments}
  \item{sigma_mu}{spread of mus when proposing new segments}
  \item{subst_hyper_mean}{prior mean of log kappa}
  \item{subst_hyper_variance}{prior variance of log kappa}
  \item{diver_hyper_mean}{prior mean of log mu}
  \item{diver_hyper_variance}{prior variance of log mu}
  \item{top_breaks}{1 = include topology break-points into the model, 0 = do not include}
  \item{par_breaks}{1 = include substitution change-points in the model, 0 = do not include}
  \item{step.size}{a positive integer indicating the number of base pairs to move the window used to create the set of candidate trees}
  \item{boot}{a non-negative integer indicating the number of times to bootstrap each window when creating candidate trees}
}
\details{DualBrothers is a recombination detection software based on the dual Multiple Change-Point model. This model allows for changes in topology and evolutionary rates across sites in a multiple sequence alignment. DualBrothers uses a Bayesian approach together with an MCMC sampling to simulate from the posterior distribution of the dual MCP model parameters. More information can be found at http://www.biomath.ucla.edu/msuchard/DualBrothers/.}
\value{DualBrothers creates a basename.post file in the current working directory.}
\references{Minin VN, Dorman KS, Suchard MA. (2005). Dual multiple change-point model leads to more accurate recombination detection. Bioinformatics, 21:3034-3042.}
