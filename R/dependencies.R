# dependencies


url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
pkgFile <- "JohnsonDistribution_0.24.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)

# source install earlier version of gls for multinma
url <- "https://cran.r-project.org/src/contrib/Archive/gsl/gsl_2.1-8.tar.gz"
pkgFile <- "gsl_2.1-8.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)


library(gridExtra)

library(tidyverse)
library(netmeta)
library(emmeans)
library(multinma)
library(twang)

#library(optweight)

library(devtools)
library(remotes)

pak::pak("bonorico/gcipdr")

library(gcipdr)