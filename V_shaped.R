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
                           levels = c("none", "mild", "high"))
  )


# plot 

ggplot(
  res1,
  aes(x = samplesize, y = estimate, colour = contrast, shape = evidence)
) +
  geom_point(size = 3) +
facet_wrap(vars(inconsistency)) +
  geom_point()
  xlab("Sample size") +
  ggtitle("")

