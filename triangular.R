# simulate triangular network ABC and corresponding ideal multiarm trial (IMT)
# 

source("./R/trial_simul2.R")
library(tidyverse)
library(netmeta)

#### IMT ####
# settings:
# N = 10000 large enough ideal trial
# average TRT effects = BA: -10; CA = -5 
# modifier = two balanced strata P(V1)=0.5, P(V2)=0.5, 
# treatment effect in V1 = 0
# treatment effect in V2 = BA: -20; CA: -10
# indirect effect BC, by consistency = BA-CA= -10 - (-5) = -5

imt <- trial_simul2(
  N = 10000,
  delta = c(-10, -5),
  mod_dist = 0.5,
  deltasub = 0
)

# check

imtmod <- lm(y~trt_name + x, data = imt$data)

summary(imtmod)
# check consistency
update(imtmod, contrasts = list(
  trt_name = contr.treatment(
    n=LETTERS[1:3], base = 3)))

         
## simulate network
## settings:
## N studies = 5
## designs = BA, CA (head-to-head only, no multiarm)
## each study has same the prognostic and modifier distribution as IMT
## common effect = all effects are assumed fixed across different study populations
## sigma = residual variance can vary from study to study
## 


set.seed(45)
# network settings
settings <- list(

  N = round(runif(5, min = 100, max = 500)),
  design = list(
    c("A", "B"), c("A", "B"), c("A", "B"),
    c("A", "C"), c("A", "C")
  ),
  delta = c(rep(-10, 3), rep(-5, 2)),
  subdelta = rep(0, 5),
  mod_prev = rep(0.5, 5),
  sigma = runif(5, 0.5, 3) # study specific precision
    
)


# stack network data

ipd_network <- lapply(1:5,
       function(i)
         trial_simul2(
           N = settings$N[i],
           trt_names = settings$design[[i]],
           delta = settings$delta[i],
           deltasub = settings$subdelta[i],
           mod_dist = settings$mod_prev[i],
           sigma0 = settings$sigma[i],
           seed = 12+i
         )$data |> 
         add_column(study = as.character(i),
                    .before = 1)
) |> 
  bind_rows()


# estimand: average treatment effect. It is assumed that there is no prior knowledge of the effect modification (otherwise it would be most likely controlled for)

# estimate: calculate average treatment effect for each study by adjusting for known prognostic factor BUT not for effect modifier
split(ipd_network, ipd_network$study) |>
  map_df(
    \(x) {
      mod <- lm(y~trt_name + x, data = x)
      
    } 
  )
  
