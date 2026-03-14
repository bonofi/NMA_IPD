# simulate triangular network ABC and corresponding ideal multiarm trial (IMT)
# 

source("./R/trial_simul2.R")

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

         
