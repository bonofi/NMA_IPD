#' Balance studies using IPW
#' 
#' 


ipw_balance <- function(ipd_network, 
                        model_formula = as.formula(study ~ V1 + V2),  # do not balance for X
                        estimand = c("ATE", "ATT"),
                        stop_rule = "ks.mean",   # can be a vector
                        n_trees = 3000)
{
  
  estimand <- match.arg(estimand)
  
  browser()
  
  # factorize study label
  ipd_network$study <- as.factor(ipd_network$study)

  
  res <- twang::mnps(
    model_formula,
    data = ipd_network,
    estimand = estimand,
    verbose = FALSE,
    stop.method = stop_rule,
    n.trees = n_trees)
 
  
  # check convergence
  #par(mfrow=c(3, 2))
 
   res$psList[[1]]$ps[res$psList[[1]]$treat == 1, 1] |> hist()
 
   res$psList[[1]]$w[res$psList[[1]]$treat == 1, 1] |> hist()
   
   bal.table(res) |> 
     filter(var == "V1")
   
   plot(res$psList[[1]])
}


ipw_balance(
  res1dat |> 
    filter(
      inconsistency == "high",
      samplesize == "small"
    )
)