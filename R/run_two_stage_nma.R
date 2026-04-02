#' Run Two-stage NMA
#' 
#' 



run_two_stage_nma <- function(
    ipd_network,  # coming from trial_simul2 
    study_level_model_formula = formula(y~trt_name + x),
    nma_common = TRUE
    ){
  
  
  # Stage 1 
  raw_lm <- split(ipd_network, ipd_network$study) |>
    map_df(
      \(df) lm(
        study_level_model_formula, 
        data = df) |> 
        broom::tidy() |> 
        add_column(treat2 = df$ref_trt |> 
                     unique() |> 
                     stringr::str_remove_all("seq:0; name:"),
                   studlab = unique(df$study))
    )
  
  
  netdata <- raw_lm |> 
    filter(grepl("trt_", term)) |> 
    mutate(
      term = stringr::str_remove_all(term,
                                     "trt_name")
    ) |> 
    select(term, estimate, std.error, treat2, studlab) |> 
    rename(
      treat1 = term,
      TE = estimate,
      seTE = std.error
    )
  
  # Stage 2
  nma <- netmeta(TE, seTE, treat1, treat2, 
                 studlab, data = netdata)
  
  netmeta::netgraph(nma)
  
  # inconsistency test
  netmeta:::forest.netsplit(
    netmeta::netsplit(nma, show = "all")
  )
  
  nmatype <- ifelse(nma_common, "common", "random")
  message(paste0(
    "Two-stage NMA: ", nmatype, " effect estimates"
  ))
  # to invert sign use lower.tri
  nmacontr <- data.frame(
    contrast = combn(
      rownames(nma[[paste0("TE.", nmatype)]]), 2) |> 
      apply(2, paste, collapse = " - "),
    estimate = nma[[paste0(
      "TE.", nmatype)]][upper.tri(nma[[paste0("TE.", nmatype)]])],
    SE = nma[[paste0("seTE.", nmatype)]][upper.tri(nma[[paste0("TE.", nmatype)]])],
    df = NA,
    stat = nma[[paste0("statistic.", nmatype)]][upper.tri(nma[[paste0("TE.", nmatype)]])],
    p.value = nma[[paste0("pval.", nmatype)]][upper.tri(nma[[paste0("TE.", nmatype)]])],
    lower.CL = nma[[paste0("lower.", nmatype)]][upper.tri(nma[[paste0("TE.", nmatype)]])],
    upper.CL = nma[[paste0("upper.", nmatype)]][upper.tri(nma[[paste0("TE.", nmatype)]])],
    model = paste0(nmatype, " effect")
  )
  
  cat(
    " Number of studies: k =", nma$k, "\n",
    "Number of pairwise comparisons: m =", nma$m, "\n",
    "Number of treatments: n =", nma$n, "\n",
    "Number of designs: d =", nma$d, "\n"
  )
  
  print(nmacontr |> 
          dplyr::select(
            !dplyr::any_of("model")
          ) |>  
          dplyr::mutate(
            dplyr::across(
              dplyr::where(is.numeric), 
              \(x) round(x, 2)
            )
          )
  )  
  
  
  return(
    list(
      NMA = nma,
      table = nmacontr |> 
        tibble::add_column(
        evidence = "2-stage adjusted NMA", #adjusting for moderator variable
        estimand = "ATE"
      )
    )
  )
  
}