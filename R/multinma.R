#' MultiNMA balancing
#' 
#' 
#' 


multinma <- function(ipd_network, 
                     modelformula = as.formula(~x + V),   # ~(x + V)*.trt
                     datalevel = c("ipd", "agd", "ipd-agd"),
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

  # must change name of prognostic variable due to bugs
  
  ipd_network <- ipd_network |> 
    dplyr::rename(X = x)
  
  replform <- modelformula[2] |> as.character() |> 
    stringr::str_replace("x", "X")
  
  modelformula <- paste0("~", replform)
  
  ### prepare AGD in case needed
  
  agd_network <- ipd_network |> 
    dplyr::select(study, trt_name, y, X, 
                  starts_with("V")) |> 
    dplyr::group_by(study, trt_name) |> 
    dplyr::summarise_if(
      is.numeric,
      list(mean=mean, sd = sd), na.rm = TRUE
    )  |> 
    dplyr::left_join(
      ipd_network |> 
        dplyr::group_by(study, trt_name) |> 
        dplyr::summarise(n = n()),
      by = c("study", "trt_name")
    ) |> 
    dplyr::ungroup()
  
  
  if (datalevel == "ipd")
    
    pso_net <- combine_network(
      set_ipd(ipd_network, 
              study = study, 
              trt = trt_name, 
              y = y,
              trt_ref = "A"
      )
    )
  
  else if (
    datalevel == "ipd-agd"
  ) {
    # put all studies except study 1 in summary arm-level format
    # combine AGD and IPD
    
    
    pso_net <- combine_network(
      set_ipd(
        ipd_network |> 
          dplyr::filter(study == "1"),
        study = study, 
        trt = trt_name, 
        y = y,
        trt_ref = "A"
      ),
      set_agd_arm(
        agd_network |> 
          dplyr::filter(study !="1"),
        study = study, 
        trt = trt_name, 
        y = y_mean,
        se = y_sd,
        trt_ref = "A",
        sample_size = n
      )
      
    )
    
    # integrate over distribution of reference study: this will not work with more strata: 
    # consider to give command manually as argument
    pso_net <- add_integration(
      pso_net,
      X= distr(qgamma, mean = X_mean, sd = X_sd),
      V2 = distr(qbern, prob = V2_mean),
      n_int = 1000
    )
    
    
  } else if (datalevel == "agd")
    
    pso_net <- combine_network(
      set_agd_arm(
        agd_network,
        study = study, 
        trt = trt_name, 
        y = y_mean,
        se = y_sd,
        trt_ref = "A",
        sample_size = n
      )
      
    )
  
  
  #######  ATT by predicting in new population ???
  #######  use newpop in newdata in marginal_effects
  
  # # use correlation matrix of reference study
  # corref <- ipd_network |> 
  #   dplyr::filter(study =="1") |> 
  #   dplyr::select(X, starts_with("V")) |>
  #   dplyr::select(where(is.numeric)) |>
  #   # drop first stratum
  #   dplyr::select(!any_of("V1")) |>  
  #   cor(method = "spearman")
  # 
  # newpop <- agd_network |> 
  #   filter(study == "1")
  # 
  # newpop <- add_integration(
  #   newpop,
  #   X= distr(qgamma, mean = X_mean, sd = X_sd),
  #   V2 = distr(qbern, prob = V2_mean),
  #   cor = corref,
  #   cor_adjust = "spearman",
  #   n_int = 1000
  # )
  # 
  
    
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
      estimand = "ATE",  # not aware ATT estimand possible in this outcome-regression context
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
  
}



# 
# # 
# prova <- multinma(
#   res1dat |>
#     filter(
#       inconsistency == "high",
#       samplesize == "small"
#     ),
#   modelformula = as.formula(~x + V2) #as.formula(~x + V:.trt)
# )




