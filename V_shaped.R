# V-shaped network


source("./R/trial_simul2.R")
source("./R/network_simul.R")
library(tidyverse)
library(netmeta)
library(emmeans)
#library(multinma)

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


ssizes <- list(small = 1,
               medium = 10,
               large = 100)

inconsistency <- list(
  none = rep(0.5, 5),
  mild = c(0.5, 0.8, 0.5, 0.5, 0.5),
  high = c(0.5, 0.8, 0.5, 0.8, 0.8) # more inconsistency in AC minority design
)

res1 <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      
      network_simul(
        network_settings = list(
          N = c(10, 50)*ssizes[[j]], 
          design = list(
            c("A", "B"), c("A", "B"), c("A", "B"),
            c("A", "C"), c("A", "C")
          ),
          delta = as.list(c(rep(-10, 5-2), rep(-5, 5-3))),
          subdelta = as.list(rep(0, 5)),
          mod_prev = as.list(inconsistency[[i]]),
          sigma = c(0.5, 3)
        )
      )$est |> 
      tibble::add_column(
        inconsistency = i,
        samplesize = j
      )
    
  )
) |> 
  dplyr::bind_rows() |> 
  dplyr::mutate(
    samplesize = factor(samplesize, 
                        levels = c("small", "medium", "large")),
    inconsistency = factor(inconsistency,
                           levels = c("none", "mild", "high")),
    estimate = ifelse(
      inconsistency != "none" & evidence == "IMT", # IMT result redundant in other panels
      NA_real_,
      estimate
      )
  )


# plot 

ggplot(
  res1,
  aes(x = samplesize, y = estimate, colour = contrast, shape = evidence)
) +
  geom_point(size = 3) +
facet_wrap(vars(inconsistency)) +
  geom_hline(yintercept = 10, colour  ="red", linetype = 2, alpha = 0.5) + 
  geom_hline(yintercept = 5, colour  ="red", linetype = 2)  +
  geom_hline(yintercept = -5, colour  ="red", linetype = 2) +
  xlab("Sample size") +
  ggtitle("Estimation in a V-shaped treatment network for increasing inconsistency levels (Nr of studies: 5)")


# INTERPRETATION: 
# 1) inconst none. In medium sample size, bias in AB contrast is entirely due to random error (random precision and/or residual error). This alone is sufficient to cause inconsistency in BC contrast. As precision becomes uniformly high (large sample size), NMA estimates are unbiased.
# 2) inconsistent. Weight of direct/indirect evidence, which design is affected by bias and how many studies, all determine how bias propagates across the network.Focusing on large sample size, inconsistency affects only one study in majority design AB (mild scenario) also causing bias in BC contrast. Interestingly, by increasing bias for both AB and AC designs seems to cancel out bias in BC design. However, this should be regarded as a change outcome and it stresses how unpredictable bias propagation can be as bias and network size increases.  