#' Balance studies using IPW
#' 
#' 


ipw_balance <- function(ipd_network, 
                        model_formula = as.formula(study ~ x + V1 + V2),  # do not balance for X
                        estimand = c("ATE", "ATT"),
                        stop_rule = "ks.mean",   # can be a vector
                        n_trees = 3000)
{
  
  estimand <- match.arg(estimand)

  
  # factorize study label
  ipd_network$study <- as.factor(ipd_network$study)

  
  res <- twang::mnps(
    model_formula,
    data = ipd_network,
    estimand = estimand,
    verbose = FALSE,
    stop.method = stop_rule,
    n.trees = n_trees)
  
  browser()
  
  # res$psList[[1]]$ps[res$psList[[1]]$treat == 1, 1] |> hist()
  # 
  # res$psList[[1]]$w[res$psList[[1]]$treat == 1, 1] |> hist()
  # 
  # bal.table(res) |> 
  #   filter(var == "V2")

  n <- length(res$psList)
  weights <- lapply(1:n,
                    function(i)
                      data.frame(
                        ps = res$psList[[i]]$ps[res$psList[[i]]$treat == 1, 1],
                        ps_weight = res$psList[[i]]$w[res$psList[[i]]$treat == 1, 1],
                        study = i
                      )
  )|> 
    dplyr::bind_rows() |> 
    tibble::rownames_to_column("subjid") |> 
    dplyr::mutate(study = as.factor(study))
  
  
  ggplot(weights,
         aes(study, ps)) + 
    geom_boxplot()
  
  data <- ipd_network |> 
    # not super clean subject identification but should work (optimal i-studyid)
    tibble::rownames_to_column("subjid") |> 
    dplyr::left_join(
      weights,
      by=c("study", "subjid")
    )
  
  # run model with weights
  
  mod <- lm(y~trt_name, data = data, weights = ps_weight)
  
}


### Plot tools for diangostics

mnps_plot <- function(res, # mnps object
                      layer = 1){
  
  n <- length(res$psList)
  k <- round(n/2)
  
  diagnostics <- c(
    "1" = "Balance convergence: ",
    "2" = "PS overlap: ",
    "3" = "Before-After weighting: ",
    "4" = "Before-After weighting: ",
    "5" = "Before-After weighting: "
  )
  

  
  # convergence
  par(mfrow = c(k, k))
  for (i in 1:n)
    plot(
      res$psList[[i]],
      plots = layer,
      main  =  paste(ifelse(i == 1, 
                            diagnostics[i], 
                            ""), 
                     "Study ", i, "vs others")
    )
  
  
}

prova <- ipw_balance(
  res1dat |> 
    filter(
      inconsistency == "high",
      samplesize == "small"
    )
)