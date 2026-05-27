#' check GC performance with IPW

ref <- res1dat |> filter(inconsistency == "high" & samplesize == "large")
dat <- rawGC4$high$large$pseud[[1]]

dat <- rawGCtest2b$pseud[[5]]

provaref <- ipw_balance(
  ipd_network = ref,
  model_formula = as.formula(study ~ x + V1 + V2),
  estimand = "ATT",
  stop_rule = "es.mean"
)

# 
prova <- ipw_balance(
  ipd_network = dat ,#|> filter(study != "5"),
  model_formula = as.formula(study ~ x + V1 + V2),
  estimand = "ATT",
  stop_rule = "ks.mean",
  n_trees = 3000
)

prova$rawres |> print()

# effect in study 3, while consistent with 1, is not properly reproduced in GC ... But after discarding 3, it appears to GC cannot well balance across studies --> residual inconsistency
prova2 <- ipw_balance(
  ipd_network = dat |> filter(study %in% c("1", "5", "4")),
  model_formula = as.formula(study ~ x + V1 + V2),
  estimand = "ATT",
  stop_rule = "ks.mean"
)







####### rerun GC with different V level 

system.time(
  rawGCtest3 <- do_gcipdr(
    ipd_network = ref,
    boot_iter = 10,
    method = "4",
    SI_k = 60000,
    only_SI = TRUE,
    seed = 30697,
    cores = 6,
    drop_ref_V = "V2"
  )
)
#########




####### try to model trt*V interaction in GC --> WRONG!!!

system.time(
  rawGCtest2b <- do_gcipdr(
    ipd_network = ref |> 
      rowwise() |> 
      # must model with complementary level of V to avoid collinearity
      mutate(inter = trt*V1) |> 
      ungroup(),
    boot_iter = 10,
    method = "3",
    SI_k = 60000,
    only_SI = TRUE,
    seed = 30697,
    cores = 6
  )
)
#########


####### try to model trt*V interaction in GC avoid redundancy

system.time(
  rawGCtest2b <- do_gcipdr(
    ipd_network = ref |> 
      rowwise() |> 
      # must model with complementary level of V to avoid collinearity
      mutate(inter = trt*V1) |> 
      ungroup(),
    boot_iter = 10,
    method = "3",
    SI_k = 60000,
    only_SI = TRUE,
    seed = 30697,
    cores = 6,
    interaction_only = TRUE
  )
)


rawGCtest2b$pseud <- lapply(
  rawGCtest2b$pseud,
  function(x) 
    x |> 
    rowwise() |> 
    mutate(
      trt = case_when(
        inter == 1 ~ 1,
        inter == 0 ~ rbinom(1,1,0.5)
      ),
      V1 = ifelse(inter == 1, 1, 1-trt),
      V2 = 1-V1
    ) |>
    ungroup() |> 
    left_join(
      ref |> 
        distinct(study, trt, trt_name),
      by = c("study", "trt")
    )
)

#########




res <- twang::mnps(
  as.formula(study ~ x + V1 + V2),
  data = dat |> mutate(study = as.factor(study)),
  estimand = "ATT",
  treatATT = "1",
  interaction.depth = 10,
  verbose = FALSE,
  stop.method = "es.mean",
  n.trees = 3000)

## weighted regression according to twang documentation: SE are different
dat$w <- twang::get.weights(res, stop.method = "es.mean")
datdes <- survey::svydesign(ids=~1, weights=~w, data=dat)

glm1 <- survey::svyglm(y ~ as.factor(trt_name), design = datdes)
summary(glm1)

