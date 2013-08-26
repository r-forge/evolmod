\name{conv.evol}
\alias{conv.evol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{blah}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
conv.evol(alignment, perc, taxa1, taxa2, map)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{
%%     ~~Describe \code{alignment} here~~
}
  \item{perc}{
%%     ~~Describe \code{perc} here~~
}
  \item{taxa1}{
%%     ~~Describe \code{taxa1} here~~
}
  \item{taxa2}{
%%     ~~Describe \code{taxa2} here~~
}
  \item{map}{
%%     ~~Describe \code{map} here~~
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
function (alignment, perc, taxa1, taxa2, map) 
{
    stop.codons <- c("TAA", "TAG", "TGA")
    nucl <- c("A", "G", "C", "T")
    align.taxa.conv <- alignment
    taxa1.aa <- codon.map[match(alignment[taxa1, ], map[, 1]), 
        2]
    taxa2.aa <- codon.map[match(alignment[taxa2, ], map[, 1]), 
        2]
    ind <- which(taxa1.aa != taxa2.aa)
    num.diff <- dim(alignment[c(taxa1, taxa2), ind])[2]
    for (i in 1:num.diff) {
        aa.1 <- unlist(strsplit(alignment[taxa1, ind[i]], ""))
        aa.2 <- unlist(strsplit(alignment[taxa2, ind[i]], ""))
        diff <- which(aa.1 != aa.2)
        if (length(diff) == 1) {
            if (rbinom(1, 1, perc) == 1) {
                free <- nucl[nucl != aa.1[diff] & nucl != aa.2[diff]]
                test.split1 <- aa.1
                test.split1[diff] <- free[1]
                test.split2 <- aa.1
                test.split2[diff] <- free[2]
                test.comb1 <- paste(test.split1[1], test.split1[2], 
                  test.split1[3], sep = "")
                test.comb2 <- paste(test.split2[1], test.split2[2], 
                  test.split2[3], sep = "")
                both <- c(test.comb1, test.comb2)
                nonstop <- both[which(!is.element(both, stop.codons))]
                nonstop <- nonstop[which(!is.element(nonstop, 
                  alignment[-c(taxa1, taxa2), ind[i]]))]
                if (length(nonstop) >= 1) {
                  replace <- sample(nonstop, 1)
                  align.taxa.conv[taxa1, ind[i]] <- replace
                  align.taxa.conv[taxa2, ind[i]] <- replace
                  others <- align.taxa.conv[-c(taxa1, taxa2), 
                    ind[i]]
                  distinct <- length(unique(others))
                  if (distinct > 1) {
                    if (distinct > 2) {
                      align.taxa.conv[-c(taxa1, taxa2), ind[i]] <- sample(others, 
                        1)
                    }
                    else {
                      major <- which(colSums(sapply(others, "==", 
                        align.taxa.conv[, ind[i]])) == 2)
                      align.taxa.conv[-c(taxa1, taxa2), ind[i]] <- others[major][1]
                    }
                  }
                }
            }
        }
    }
    return(align.taxa.conv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
