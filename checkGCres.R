#' check GC performance with IPW


dat <- rawGC4$high$large$pseud[[1]]


prova <- ipw_balance(
  ipd_network = dat,
  model_formula = as.formula(study ~ x + V1 + V2),
  estimand = "ATT",
  stop_rule = "es.mean"
)

prova$rawres |> print()


prova2 <- ipw_balance(
  ipd_network = dat,
  model_formula = as.formula(study ~ x + V1 + V2),
  estimand = "ATT",
  stop_rule = "es.mean",
  n_trees = 10000
)
