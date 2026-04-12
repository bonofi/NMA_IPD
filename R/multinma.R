#' MultiNMA balancing
#' 
#' 
#' 


multinma <- function(ipd_network, 
                     modelformula = as.formula(~x + V),   # ~(x + V)*.trt
                     datalevel = c("ipd", "agd", "ipd-agd"),
                     estimand = c("ATE", "ATT"),
                     model = c("fixed", "random"),
                     reverse_effects = TRUE,   # 
                     n_chains = 4,
                     n_iter = 2000,
                     seed = 87632
                     )
  {
  
  datalevel <- match.arg(datalevel)
  model <- match.arg(model)
  estimand <- match.arg(estimand)
  
  browser()
  
  ### prepare AGD in case needed
  
  agd_network <- ipd_network |> 
    dplyr::group_by(study, trt_name) |> 
    dplyr::summarise(
      y_mean = mean(y, na.rm = TRUE),
      y_sd = sd(y, na.rm = TRUE),
      x_mean = mean(X, na.rm = TRUE),
      x_sd = sd(X, na.rm = TRUE),
      V2_mean = mean(V2, na.rm = TRUE),
      n = n()
    ) |> 
    dplyr::ungroup()
  
  
  if (datalevel == "ipd" & estimand == "ATE")
    
    pso_net <- combine_network(
      set_ipd(ipd_network, 
              study = study, 
              trt = trt_name, 
              y = y,
              trt_ref = "A"
              )
    )
  else if (datalevel == "ipd" & estimand == "ATT") # Ref study 1
  {
    pso_net <- combine_network(
      set_ipd(
        ipd_network |> 
          dplyr::filter(study != "1"),
        study = study, 
        trt = trt_name, 
        y = y,
        trt_ref = "A"
      ),
      set_agd_arm(
        agd_network |> 
          dplyr::filter(study =="1"),
        study = study, 
        trt = trt_name, 
        y = y_mean,
        se = y_sd,
        trt_ref = "A",
        sample_size = n
      )
      
    )
    
    # use correlation matrix of reference study
    corref <- ipd_network |> 
      dplyr::filter(study =="1") |> 
      dplyr::select(X, V2) |> 
      cor(method = "spearman")
    
    # integrate over covariate distribution of reference study
    pso_net <- add_integration(
      pso_net,
      X= distr(qgamma, mean = x_mean, sd = x_sd),
      V2 = distr(qbern, prob = V2_mean),
      cor = corref,
      cor_adjust = "spearman",
      n_int = 1000
    )
    
  }
  
  else if (
    datalevel == "ipd-agd"
  ) {
    # put all studies except study 1 in summary arm-level format
    # combine AGD and IPD
    
  } else if (datalevel == "agd")
    
    
  png("./output/network.png")
  multinma:::plot.nma_data(
    pso_net, weight_nodes = TRUE, 
    weight_edges = TRUE, show_trt_class = FALSE) + 
    ggplot2::theme(legend.position = "bottom", 
                   legend.box = "vertical") 
  
  dev.off()
  
  

  
  set.seed(seed)
  
  nma <- nma(pso_net, 
             trt_effects = model,
             link = "identity", 
             likelihood = "normal",
             regression = modelformula,
             iter = n_iter,
             chains = n_chains
             )
  
  
  modcontr <- multinma:::relative_effects(nma, all_contrasts = TRUE)$summary |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      contrast = ifelse(reverse_effects,
                        paste0(.trta, " - ", .trtb),
                        paste0(.trtb, " - ", .trta)
      ),
      estimate = ifelse(reverse_effects,
                        -mean,
                        mean
      ),
      lower.CL = ifelse(reverse_effects,
                        -`97.5%`,
                        `2.5%`
      ),
      upper.CL = ifelse(reverse_effects,
                        -`2.5%`,
                        `97.5%`
      ),
      df = NA,
      stat = NA,
      p.value = NA,
      model = ifelse(model == "fixed", "common effect", 
                     ifelse(model =="random", "random effect", NA)),
      evidence = "ML-NMR",  # multilevel network metaregression
      estimand = estimand,
      level = toupper(datalevel),
      evidence2 = "Balanced" 
    ) |> 
    dplyr::rename(
      SE = sd
    ) |> 
    dplyr::select(
      contrast, estimate, SE, df, stat, p.value,
      lower.CL, upper.CL, model, evidence, estimand,
      level, evidence2
    )
  
  
  return(
    list(
      mlnmr = nma,
      est = modcontr
    )
  )
  # TODO: 
  # output for IPD-AGD model
  # only AGD...
  
}



# 
# 
prova <- multinma(
  res1dat |>
    filter(
      inconsistency == "high",
      samplesize == "small"
    ) |> 
    rename(X = x),
  modelformula = as.formula(~X + V1) #as.formula(~x + V:.trt)
)



multinma:::marginal_effects(
  prova$mlnmr, 
  baseline = "1", 
  newdata = data.frame(x = 60, V = rep(c("level1", "level2"), c(100,100))))


multinma:::predict.stan_nma(
  prova$mlnmr, 
  baseline = "1", 
  newdata = res1dat |> 
    filter(study == "1",
           inconsistency == "none",
           samplesize == "large"
           ) |>
    select(x, V),
  level = "individual"
    )


ref1 <- res1dat |> 
  filter(study == "1",
         inconsistency == "none",
         samplesize == "large"
  ) |>
  rename(X = x) |> 
  select(X, V1) |> 
  summarise(
    x_mean = mean(X),
    x_sd = sd(X),
    V1_mean = mean(V1) 
  ) |> as.data.frame()

corref <- res1dat |> 
  filter(study == "1",
         inconsistency == "none",
         samplesize == "large"
  ) |> rename(X = x) |> 
  select(X, V1) |> cor(method = "spearman")


ref1b <- add_integration(
  ref1,
  X= distr(qgamma, mean = x_mean, sd = x_sd),
  V1 = distr(qbern, prob = V1_mean),
  cor = corref,
  cor_adjust = "spearman",
  n_int = 1000
)

multinma::marginal_effects(prova$mlnmr,
                           baseline = "1",
                           newdata = ref1b)


####  ATT analysis with ML-NMR (all IPD available)




