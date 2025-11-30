#' Potential outcomes
#'
#'
#'

po_simul <-
function(m1, m0, s1 = 1, s0 = 1, pt = 0.5, N = 10, seed = 0614){

  # true ATE

  ate <- m1 - m0

  set.seed(seed)

  # define a potential outcome distribution for the untreated

  y1 <- rnorm(N, m1, s1)  # potential outcome if treated

  # define a potential outcome distribution for the treated

  y0 <- rnorm(N, m0, s0) # potential outcome if untreated

  # treatment assignment

  a <- rbinom(N, 1, pt)

  # consistency equation

  y <- y1*a + y0*(1-a)

  return(

    data.frame(
      Y1 = y1,
      Y0 = y0,
      A = a,
      delta = y1 - y0,   # individual causal effect
      Y = y
    )
  )

}
