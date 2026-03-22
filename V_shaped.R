# V-shaped network


source("./R/trial_simul2.R")
source("./R/network_simul.R")
library(tidyverse)
library(netmeta)

# GENERAL SETTINGS
##### IMT ####
# settings:
# N = 10000 large enough ideal trial
# average TRT effects = BA: -10; CA = -5 
# modifier = two balanced strata P(V1)=0.5, P(V2)=0.5, 
# treatment effect in V1 = 0
# treatment effect in V2 = BA: -20; CA: -10
# indirect effect BC, by consistency = BA-CA= -10 - (-5) = -5
## simulate network ######
## settings:
## N studies = 5
## designs = BA, CA (head-to-head only, no multiarm)
## each study has same prognostic and modifier distribution as IMT
## common effect = all effects are assumed fixed across different study populations
## sigma = residual variance can vary from study to study


res1 <- network_simul()

