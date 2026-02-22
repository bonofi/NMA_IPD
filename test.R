





po <- po_simul(m1 = 10, m0 = 20)

# true ate
10 - 20

mean(po$delta)

sd(po$delta)


mean(
  po$Y[po$A == 1]

) - mean(
  po$Y[po$A == 0]
)



############### INTERACTION #################


interaction <- function(n = 1000,
                        beta0 = 1,
                        beta1 = 2,
                        betatrt = -10,
                        betaint = 5,
                        sigma0 = 1,
                        fx = function(x) rgamma(x, 60) # rbinom(x, 1, 0.)
                        )
{
  
  set.seed(12)
  
  n <-  1000
  
  # simulate predictor
  x <- fx(n)
  
  # random allocation
  a <- rbinom(n, 1, 0.5)
  
  # mu
  mu <- beta0 + beta1*x + betaint*x*a + betatrt*a
  
  y <- rnorm(n, mu, sigma0)
  
  data.frame(y=y, trt=a, x=x)
  
}


dd <- interaction()

# change predictor dist
dd2 <- interaction(
  fx = function(x) rgamma(x, 20)
)


summary(lm(y~trt, data = dd))

summary(lm(y~trt + x, data = dd))

summary(lm(y~x*trt, data = dd))


