#' Simulate Trial
#' @description
#' Simulate multi arm (arms >=2) trial from IPD generating mechanism
#' @description
#' Simulate trial with continuous primary outcome and residual unit-variance Normal error, e.g. cholesterol levels, in presence of prognostic variable and its interaction with the treatment effect.
#'
#' @param pt - allocation ratio, 1:1 -> 0.5
#' @param fx - function to simulate prognostic variable vector, e.g. age with mean 60, Gamma(60)
#' @param deltasub - vector of treatment effect for K-1 ordered subgroups
#' @param mod_dist - vector of prevalence of K subgroups (categorical effect modifier)


trial_simul2 <- function(N, delta, mu0, beta, deltasub = c(0, 10), sigma0 = 1, pt = 0.5,
                        fx = function(x) rgamma(x, 60), mod_dist = c(0.2, 0.4, 0.4),
                        seed = 5602783, trt_names = LETTERS[1:(length(delta) + 1)]){
  
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
  
  browser()
  #  mean of deltas equal delta (total trt effect)
  # deltas <- c(
  #   deltasub,
  #   (delta*dim(z)[2]) - sum(deltasub)   # derive effect of third subgroup
  # )
  
  deltas <- rbind(
    deltasub,
    matrix( 
      (delta*dim(z)[2]) - apply(deltasub, 2, sum),
      nrow  =1)
  )
  
  # actual subgroup effect based on subgroup
  subdelta <- z%*%deltas
  
  if (is.na(pt))
    #equal allocation
    pt <- (N/(length(delta) + 1))/N
  
  # generate allocation list and 
  # Outcome expectation: linear relationship with effect modification
  if ( (length(delta) == length(pt)) == 1 )
  {
    a <- rbinom(N, 1, pt)
    y <- y0 + subdelta*a
    
  } else {
    a <- rmultinom(N, 1, prob = pt) |>
      t()
    # drop reference (placebo) arm for calculating sub effects
    y <- y0 + (subdelta*a[, -1]%*%rep(1,length(delta))) 
  }
     

  # simulate one head-to-head trial
  
  trt_seq <- sort(unique(
    a%*%c(1:dim(a)[2]) - 1
  )) |> 
    as.vector()
  
  trt_label <- data.frame(
    trt = trt_seq,
    trt_name = sort(unique(trt_names))
  )
  
  
  out <- tibble::tibble(
    epsilon = res_err,
    mu0 = mu0,
    beta = beta,
    tot_delta = delta,
    sub_delta = subdelta,
    y = y,
    trt = trt_seq,
    x = x,
    V = paste(
      "level",
      z%*%c(1:dim(z)[2]),
      sep ="_"
    )
  ) |>
    cbind(
      as.data.frame(z)
    ) |> 
    dplyr::left_join(
      trt_label,
      by = "trt"
    )
  
  
  return(out)
}

# test

dat <- trial_simul2(N = 1000, delta = -10, mu0 = 20, beta = 2)

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

