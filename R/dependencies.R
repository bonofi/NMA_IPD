# dependencies


url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
pkgFile <- "JohnsonDistribution_0.24.tar.gz"
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

#library(devtools)
#library(remotes)

#devtools::install_github("bonorico/gcipdr")
#pak::pak("bonorico/gcipdr")

#library(gcipdr)

# source files manually
# sapply(
#   list.files("./gcipdr/R/", 
#              full.names = TRUE), source)