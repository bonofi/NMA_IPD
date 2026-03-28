# test trial_simul

source("./R/trial_simul.R")

dat <- trial_simul(N = 1000, delta = -10, mu0 = 20, beta = 2)

# change predictor distribution: younger people
dat2 <- trial_simul(N = 1000, delta = -10, mu0 = 20, beta = 2, 
                    fx = function(x) rgamma(x, 30))

# change effect modifier distribution: younger people
dat3 <- trial_simul(N = 1000, delta = -10, mu0 = 20, beta = 2, 
                    mod_dist = c(0.1, 0.1, 0.8))


# estimator

# average effect is slightly biased due to unequal prevalence of the modifier strata. This bias can be addressed by adjusting for the modifier (via subgroup analysis or interaction term) to obtain the unbiased stratum specific effect whose average is the true average trt effect. Any prognostic variable must be accounted for to increase precision
summary(
  lm(y~trt, data = dat)
)

# subgroup analysis returns true subgroup effects but slighly biased if prognostic factor, X, not accounted for

summary(
  lm(y~trt + x, data = dat, subset = V1 == 1)
)

summary(
  lm(y~trt + x, data = dat, subset = V2 == 1)
)

summary(
  lm(y~trt + x, data = dat, subset = V3 == 1)
)


#################

# including an interaction with the moderator and adjusting for prognostic factor X is akin to performing the subgroup analysis above ...

summary(
  lm(y~trt + V, data = dat)
)

## 
summary(
  lm(y~trt + V + x, data = dat)
)
##
summary(
  lm(y~trt + V*trt, data = dat)
)
###


# the interaction terms are the V2-V1 = 10 and V3-V1 = -40 contrasts respectively. Summing the interaction term with trt yields the stratum specific effect (trt = effect of V1). 
summary(
  lm(y~trt + V*trt + x, data = dat)
)


# decomposition
summary(
  lm(y~trt + x*trt, data = dat)
)

summary(
  lm(y~trt + x, data = dat)
)

# subgroup analysis

lm(y~trt, data = dat, subset = x < 54)$coef
lm(y~trt, data = dat, subset = x > 54 & x < 64)$coef
lm(y~trt, data = dat, subset = x > 64)$coef

