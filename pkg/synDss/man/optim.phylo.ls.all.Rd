\name{optim.phylo.ls.all}
\alias{optim.phylo.ls.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{another blah}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
optim.phylo.ls.all(D, set.neg.to.zero = TRUE, fixed = FALSE, tol = 1e-10, collapse = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{set.neg.to.zero}{
%%     ~~Describe \code{set.neg.to.zero} here~~
}
  \item{fixed}{
%%     ~~Describe \code{fixed} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{collapse}{
%%     ~~Describe \code{collapse} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (D, set.neg.to.zero = TRUE, fixed = FALSE, tol = 1e-10, 
    collapse = TRUE) 
{
    if (class(D) == "dist") 
        D <- as.matrix(D)
    n <- nrow(D)
    all.trees <- allTrees(n, tip.label = row.names(as.matrix(D)))
    allQ <- vector()
    bestQ <- Inf
    for (i in 1:length(all.trees)) {
        all.trees[[i]]$edge.length <- rep(dim(all.trees[[i]]$edge)[1], 
            1)
        all.trees[[i]] <- nnls.tree(D, all.trees[[i]], trace = 0)
        allQ[i] <- attr(all.trees[[i]], "RSS")
    }
    best <- which(allQ == min(allQ))
    if (length(best) > 1) {
        best <- best[1]
    }
    return(all.trees[[best]])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
