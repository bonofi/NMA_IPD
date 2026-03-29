#' Balance studies using IPW
#' 
#' 


ipw_balance <- function(ipd_network)
{
  browser()
  
  # factorize study label
  ipd_network$study <- as.factor(ipd_network$study)

  
  res <- twang::mnps(
    study ~ x + V1 + V2,
    data = ipd_network,
    estimand = "ATE",
    verbose = FALSE,
    stop.method = c("es.mean", "ks.mean"),
    n.trees = 3000)
 
  
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