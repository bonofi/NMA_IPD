#' Simulate Trial
#' @description
#' Simulate head-to-head trial from IPD generating mechanism
#' @description
#' Simulate trial with continuous primary outcome and residual unit-variance Normal error, e.g. cholesterol levels, in presence of prognostic variable and its interaction with the treatment effect.
#'
#' @param pt - allocation ratio, 1:1 -> 0.5
#' @param fx - function to simulate prognostic variable vector, e.g. age with mean 60, Gamma(60)
#' @param deltasub - vector of treatment effect for K-1 ordered subgroups
#' @param mod_dist - vector of prevalence of K subgroups (categorical effect modifier)


trial_simul <- function(N, delta, mu0, beta, deltasub = c(0, 10), sigma0 = 1, pt = 0.5,
                        fx = function(x) rgamma(x, 60), mod_dist = c(0.2, 0.4, 0.4),
                        seed = 5602783){

  set.seed(seed)


  res_err <- rnorm(N, sd = sigma0)

  # simulate prognostic variable
  x <- fx(N)

  # generate baseline population (target population) as the combination of average outcome mu0 and prognostic variable

  y0 <- mu0 + beta*x + res_err

  # distribution of effect modifier (non prognostic)

  z <- rmultinom(N, 1, prob = mod_dist) |>
    t()

  # split total trt effect into subgroup modified effects

  # mean of deltas equal delta (total trt effect)
  deltas <- c(
    deltasub,
    (delta*dim(z)[2]) - sum(deltasub)   # derive effect of third subgroup
  )

  # actual subgroup effect based on subgroup
  subdelta <- z%*%deltas

  # generate allocation list

  a <- rbinom(N, 1, pt)

  # simulate one head-to-head trial

  # Outcome expectation: linear relationship with effect modification

  y <- y0 + subdelta*a

  out <- tibble::tibble(
    epsilon = res_err,
    mu0 = mu0,
    beta = beta,
    tot_delta = delta,
    sub_delta = subdelta,
    y = y,
    trt = a,
    x = x,
    V = paste(
      "level",
      z%*%c(1:dim(z)[2]),
      sep ="_"
    )
  ) |>
    cbind(
      as.data.frame(z)
    )


  return(out)
}

# test

dat <- trial_simul(N = 1000, delta = -10, mu0 = 20, beta = 2)

# estimator

# total effect
summary(
  lm(y~trt, data = dat)
)

# subgroup

summary(
  lm(y~trt, data = dat, subset = V1 == 1)
)

summary(
  lm(y~trt, data = dat, subset = V2 == 1)
)

summary(
  lm(y~trt, data = dat, subset = V3 == 1)
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

