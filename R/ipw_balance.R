#' Balance studies using IPW
#' 
#' 


ipw_balance <- function(ipd_network, 
                        model_formula = as.formula(study ~ x + V2),  # do not balance for X
                        estimand = c("ATE", "ATT"),
                        stop_rule = "ks.mean",   # can be a vector
                        n_trees = 3000)
{
  
  estimand <- match.arg(estimand)

  
  # factorize study label
  ipd_network$study <- as.factor(ipd_network$study)

  
  res <- twang::mnps(
    model_formula,
    data = ipd_network |> 
      as.data.frame(),
    estimand = estimand,
    verbose = FALSE,
    stop.method = stop_rule,
    n.trees = n_trees)
 

  baltable <- twang::bal.table(res)
  
  n <- length(res$psList)
  weights <- lapply(1:n,
                    function(i)
                      data.frame(
                        ps = res$psList[[i]]$ps[res$psList[[i]]$treat == 1, 1],
                        ps_weight = res$psList[[i]]$w[res$psList[[i]]$treat == 1, 1],
                        study = i,
                        usubjid = paste0(i, "-", 1:sum(res$psList[[i]]$treat == 1))
                      )
  )|> 
    dplyr::bind_rows() |> 
    dplyr::mutate(study = as.factor(study))
  
  # print diagnostics
  
  if (!dir.exists("./output/")) 
    dir.create("./output/")
  
  pdf(paste0(
    "./output/inconsistency-", unique(ipd_network$inconsistency),
    "_samplesize-", unique(ipd_network$samplesize), ".pdf"
  ))
  
  plot(res)
  plot(res, plots = 2)
  ggplot(weights,
         aes(study, ps)) + 
    geom_boxplot() + 
    labs(title = "Propensity scores")
  plot(res, plots = 3)
  plot(res, plots = 4)
  plot(res, plots = 5)
  
  grid::grid.newpage()
  # gridExtra::grid.table(
  #   baltable |> 
  #     dplyr::filter(var == "V1")
  # )
  # grid::grid.newpage()
  gridExtra::grid.table(
    baltable |> 
      dplyr::filter(var == "V2")
  )
  dev.off()
  
  
  
  data <- ipd_network |> 
    dplyr::left_join(
      weights,
      by=c("study", "usubjid")
    )
  
  browser()
  # run model with weights
  
  mod <- lm(y~trt_name, data = data, weights = ps_weight)
  
  Vpror <- glm(V2 ~ trt_name + study, data = data, weights = ps_weight) |> 
    predict(type = "response")
  
  data.frame(pred=Vpror, V2 = data$V2, study = data$study, trt = data$trt_name) |> 
    group_by(study, trt, V2) |> 
    filter(V2 ==1 & trt != "A") |> 
    summarise( mean_V2 = mean(pred)) |> 
    ungroup() |> 
    select(!any_of("V2")) |> 
    rowwise() |> 
    mutate(
      mean_V1 = 1-mean_V2,
      weight_eff = case_when(
        trt == "B" ~ weighted.mean(
          c(-20,0), c(mean_V2, mean_V1)),
        trt == "C" ~ weighted.mean(
          c(-10,0), c(mean_V2, mean_V1)),
        TRUE ~ NA_real_
        )
      ) |> 
    group_by(trt) |> 
    summarise(ATE = mean(weight_eff))
    
    
  
  
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