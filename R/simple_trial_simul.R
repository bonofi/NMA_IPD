# simple trial simulation



simple_trial_simul <- function(N, delta, mu0, beta, sigma0 = 1, pt = 0.5,
                        fx = function(x) rgamma(x, 60), confounding_x = FALSE,
                        seed = 5602783){
  
  set.seed(seed)
  
  # homoscedastic error
  res_err <- rnorm(N, sd = sigma0)
  
  # simulate prognostic variable
  x <- fx(N)
  
  # generate baseline population (target population) as the combination of average outcome mu0 
  # and prognostic variable + error
  
  y0 <- mu0 + beta*x + res_err
  
 
  # generate allocation list based on probability = pt
  
  if (confounding_x)
  {
    # ind with x below mean have Pr(a) = 0.7 else Pr(a) = 0.2
    #pt <- ifelse(x < mean(x, na.rm = TRUE), 0.7, 0.2)
    
    beta0 <- 0.1   # Intercept
    beta1 <- -0.2   # Effect of covariate
    
    # Calculate probability using logistic function
    logit <- beta0 + beta1 * (x - mean(x, na.rm = TRUE))
    pt <- 1 / (1 + exp(-logit))  # logistic transformation

    a <- sapply(pt, 
                function(x)
                  
                  rbinom(1, 1, x),
                
                simplify = TRUE
    )
  } else
    a <- rbinom(N, 1, pt)
  
  # simulate one head-to-head trial
  
  # Outcome
  
  y <- y0 + delta*a
  
  out <- tibble::tibble(
    epsilon = res_err,
    mu0 = mu0,
    beta = beta,
    tot_delta = delta,
    y = y,
    trt = a,
    x = x,
    prob_trt = pt
  )
  
  return(out)
}

# test

dat <- simple_trial_simul(N = 1000, delta = -10, mu0 = 20, beta = 2)

dat_conf <- simple_trial_simul(N = 1000, delta = -10, 
                               mu0 = 20, beta = 0.5, confounding_x = TRUE)


summary(
  lm(y~trt, data = dat)
)

# prognostic variable increases R-square
summary(
  lm(y~trt + x, data = dat)
)



summary(
  lm(y~trt, data = dat_conf)
)

summary(
  lm(y~trt + x, data = dat_conf)
)
