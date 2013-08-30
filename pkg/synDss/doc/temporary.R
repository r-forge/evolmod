setwd("C:/Users/Peter/R-Forge/pkg/synDss/doc")
library(tools)

Sweave("tutorial.Rnw")
texi2dvi("tutorial.tex",pdf=TRUE)


