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
    
  }
    
  png("./output/network.png")
  multinma:::plot.nma_data(
    pso_net, weight_nodes = TRUE, 
    weight_edges = TRUE, show_trt_class = FALSE) + 
    ggplot2::theme(legend.position = "bottom", 
                   legend.box = "vertical") 
  
  dev.off()
  
  
  # 
  # pso_net <- combine_network(
  #   set_ipd(pso_ipd, 
  #           study = studyc, 
  #           trt = trtc, 
  #           r = pasi75,
  #           trt_class = trtclass),
  #   set_agd_arm(pso_agd, 
  #               study = studyc, 
  #               trt = trtc, 
  #               r = pasi75_r, 
  #               n = pasi75_n,
  #               trt_class = trtclass)
  # )
  # 
  
  
  set.seed(seed)
  
  nma <- nma(pso_net, 
             trt_effects = model,
             link = "identity", 
             likelihood = "normal",
             regression = modelformula,
             iter = n_iter,
             chains = n_chains
             )
  
  
  browser()
  
  multinma:::relative_effects(nma, all_contrasts = TRUE)$summary |> 
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
      estimand = ifelse(
        datalevel %in% c("ipd", "agd"), "ATE", 
        ifelse(datalevel == "ipd-agd", "ATT", NA)
      ),
      level = toupper(datalevel) 
    ) |> 
    dplyr::rename(
      SE = sd
    ) |> 
    dplyr::select(
      contrast, estimate, SE, df, stat, p.value,
      lower.CL, upper.CL, model, evidence, estimand,
      level
    )
  
  
  # TODO: check TRT interactio model
  # output for IPD-AGD model
  # only AGD...
  
}





multinma(
  res1dat |> 
    filter(
      inconsistency == "high",
      samplesize == "small"
    )
)