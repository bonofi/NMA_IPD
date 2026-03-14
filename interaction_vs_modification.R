# Example showing the difference between interaction and effect modification
# Cytel Inc


# simulate interaction effect data
interaction <- function(n = 1000,
                        beta0 = 1,   #intercept
                        beta1 = 2,   # interaction variable -- define it in fx
                        betatrt = -10, # average direct TRT effect
                        betaint = 5,   # interaction with TRT term
                        sigma0 = 1,    # residual error
                        fx = function(x) rgamma(x, 60), # sample(LETTERS[1:2], x, prob = c(0.5, 0.5), replace = TRUE)      
                        seed = 9186
)
{
  
  set.seed(seed)
  
  # simulate predictor -- interaction variable
  x <- fx(n)
  
  # random allocation
  a <- rbinom(n, 1, 0.5)
  
  # mu
  mu <- beta0 + beta1*x + betaint*x*a + betatrt*a
  
  y <- rnorm(n, mu, sigma0)
  
  data.frame(y=y, trt=a, x=x)
  
}


# simulate treatment modification effect data
# effect modification occurs in subgroup A and B
modification <- function(n = 1000,
                        beta0 = 1,   #intercept
                        betatrt_average = -10, # average TRT effect over group A and B
                        betatrt_A = 20, # TRT effect in group A (effect in group B is readily obtained from average effect)
                        prob_A = 0.5,  # population prevalence of group A 
                        sigma0 = 1,    # residual error = fixed 
                        seed = 9186
                        )
{
  
  set.seed(seed)
  
  betatrt_B <- 2*betatrt_average - betatrt_A
  
  # simulate effect modifier variable
  x <- sample(LETTERS[1:2], n, 
              prob = c(prob_A, 1-prob_A), 
              replace = TRUE)
  
  # random allocation
  a <- rbinom(n, 1, 0.5)
  
  # mu
  mu <- beta0 + ifelse(x == "A", betatrt_A, betatrt_B)*a
  
  y <- rnorm(n, mu, sigma0)
  
  data.frame(y=y, trt=a, x=x)
  
}



########## UNDERSTAND THE DIFFERENCE IN THE ANALYSIS ############

# INTERACTION
# average direct effect = -10
# effect in subgroup = ? 
# average direct effect of interactor = 2
# interaction term = 5

# IT IS NOT EASY TO CONTROL THE "TRUE" SUBGROUP EFFECT USING A STATISTICAL INTERACTION FOR THE DATA GENERATING MODEL, BECAUSE THERE IS NO INTUITIVE TRANSLATION OF THE MODEL PARAMETERS (interaction term and direct effects) INTO "TRUE" AVERAGE AND SUBGROUP EFFECTS.

# continuous interaction term
dat_int <- interaction()

# average TRT effect strongly biased if interaction not accounted for !!!
summary(lm(y~trt, data = dat_int))
summary(lm(y~trt + x, data = dat_int))

# only accounting for interaction term returns unbiased trt estimates
summary(lm(y~x*trt, data = dat_int))

# SUBGROUP analysis: detects effect modification but estimates are always biased !
# on continuous scale: not clear what is subgroup cut-off to observe modification

# cutting off on mean value of interaction variable: no adjustment for interaction variable
summary(lm(y~trt, data = dat_int, subset = x < 60))
summary(lm(y~trt, data = dat_int, subset = x >= 60))
# cutting off on mean value of interaction variable: with adjustment for interaction variable
summary(lm(y~trt*x, data = dat_int, subset = x < 60))
summary(lm(y~trt*x, data = dat_int, subset = x >= 60))


## NOW using binary interaction term
dat_int2 <- interaction(fx = function(x) rbinom(x, 1, 0.5)) 
dat_int2$x <- as.factor(dat_int2$x)

# average TRT effect biased but less so on 0-1 binary scale, 
# if interaction not accounted for !!!
summary(lm(y~trt, data = dat_int2))
summary(lm(y~trt + x, data = dat_int2))

# AGAIN, only accounting for interaction term returns unbiased trt estimates
summary(lm(y~x*trt, data = dat_int2))

# SUBGROUP analysis: detects effect modification but estimates are biased --> on binary scale, the average of the group-specific effects does not return the true average effect
# on binary scale: straightforward subgroup analysis

# cutting off on mean value of interaction variable: no adjustment for interaction variable.
# Adjusting for binary interaction yields the same result
summary(lm(y~trt, data = dat_int2, subset = x == "0"))
summary(lm(y~trt, data = dat_int2, subset = x == "1"))




##### EFFECT MODIFICATION #######
# average effect = -10
# effect in A = 20
# effect in B = -40


dat_mod <- modification()

# unbiased average TRT effect 
summary(lm(y~trt, data = dat_mod))

# slightly biased average TRT effect for higher R^2 if adjusting for modifier
summary(lm(y~trt + x, data = dat_mod))

# statistical interaction highlights modification effect but returns "messed-up" estimates.
# In fact, in this clean scenario (no other prognostic variables) the TRT effect is equal to the effect in group A and the interaction is equal to the contrast between group B and A (B-A), which means the actual (unbiased) effect in group B is trt + interaction. Therefore, all parameters are thrown all over but still interpretable if put in correct slot. If prognostic factors are present and unaccounted for, the interaction terms will be biased
summary(lm(y~x*trt, data = dat_mod))

# SUBGROUP analysis: unbiased subgroup effects whose average correctly returns the true average TRT effect

summary(lm(y~trt, data = dat_mod, subset = x == "A"))
summary(lm(y~trt, data = dat_mod, subset = x == "B"))


