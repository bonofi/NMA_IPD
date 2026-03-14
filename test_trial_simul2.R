# test trial_simul2

source("./R/trial_simul2.R")

# test

dat <- trial_simul2(N = 1000, delta = -10, mu0 = 20, beta = 2)$data

# test multiarm
dat2 <- trial_simul2(N = 1000, delta = c(-10, -20), mu0 = 20, beta = 2)$data

dat2 |> dplyr::count(trt_name, V)

# test asymptotic: N -> big -> unbiased effect modifications. Average effect = weighted average of modifications with weight equal to specific distribution of modifier strata.     

dat3 <- trial_simul2(N = 10000, delta = c(-10, -20), mu0 = 20, beta = 2)$data

# ....balanced modifier strata
dat4 <- trial_simul2(N = 10000, delta = c(-10, -20), mu0 = 20, beta = 2,
                     mod_dist = c(0.333, 0.333))$data

# ....no modifier
dat5 <- trial_simul2(N = 1000, delta = c(-10, -20), mu0 = 20, beta = 2,
                     mod_dist = c(0.333, 0.333), no_modifier = TRUE)$data


# estimator

# total effect
summary(
  lm(y~trt_name, data = dat)
)

# test subgroup

summary(
  lm(y~trt_name, data = dat, subset = V1 == 1)
)

summary(
  lm(y~trt_name, data = dat, subset = V2 == 1)
)

summary(
  lm(y~trt_name, data = dat, subset = V3 == 1)
)

# test indirect effect

summary(
  lm(y~trt_name, data = dat2, subset = trt_name != "A")
)



### distribution of effect modifiers

with(dat2,
     prop.table(
       table(trt, V)
     ) 
)
# without "proportionate" stratified treatment allocation the observed frequency of V is the prevalence of V in the population, P(V), re-scaled by the allocation probability, P(A), that is, 
# P(V)*P(A). However, few trial designs have a proportionate stratification. The goal of stratified randomization is not to maintain P(V) in the trial, but to balance V across arms.

with(dat2,
     prop.table(
       table(V)
     ) 
)*0.333
