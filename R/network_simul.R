#' Simulate treatments network and corresponding ideal multiarm trial (IMT)
#' @param N - integer - size of IMT
#' @param delta - vector - TRT effects versus common reference in IMT. TRT names are automatically labelled with letters, A, B, C ...
#' @param mod_dist - vector - prevalence in TRT modifier strata in IMT
#' @param deltasub - vector - equal length of mod_dist; TRT effect in modifier's stratum in IMT  
#' @param network_settings - list - with elements: 
#' - K (number of studies) 
#' - N (vector - min-max range sample size in study; will be Uniformly sampled)
#' - design (list - TRTs in each study (head-to-head or multiarm); first name is control arm)
#' - delta (list - TRT effect(s) in each study)
#' - subdelta (list - TRT effect(s) in study's modifier stratum)
#' - mode_prev (list - prevalences of modifier's strata in each study; must have equal lenght of corresponding list element in subdelta)
#' - sigma (vector - min-max range residual error in study; will be Uniformly sampled)
#' @description
#' Number of unique TRTs in network must match number of treatments in IMT, anchored to same control. Design must be any combination of the above treatment (head-to-head or multi-arm). Design in each study determines the network structure (V-shape, triangle, star, etc ...)
#' 


network_simul <- function(
    N = 10000,
    delta = c(-10, -5),
    mod_dist = 0.5,
    deltasub = 0,
    K = 5,
    # more arguments from trial_simul2 possible (and in network setting)
    network_settings = list(
      N = c(100, 500), # by increasing N estim become more precise
      design = list(
        c("A", "B"), c("A", "B"), c("A", "B"),
        c("A", "C"), c("A", "C")
      ),
      delta = as.list(c(rep(-10, K-2), rep(-5, K-3))),
      subdelta = as.list(rep(0, K)),
      mod_prev = as.list(rep(0.5, K)),
      sigma = c(0.5, 3) # study specific precision min-max range
    ),
    rseed = 45,
    nma_common = TRUE # show NMA common or random effect results
){
  

  #### IMT ####

  imt <- trial_simul2(
    N = N,
    delta = delta,
    mod_dist = mod_dist,
    deltasub = deltasub,
    seed = rseed
  )
  
  # check
  
  imtmod <- lm(y~trt_name + x, data = imt$data)
  
  
  ## simulate network

  set.seed(rseed)
  
  Ns <- round(runif(K, 
                    min = network_settings$N[1], 
                    max = network_settings$N[2]))
  sigmas <- runif(K,
                  min = network_settings$sigma[1], 
                  max = network_settings$sigma[2])
  
  # stack network data
  
  ipd_network <- lapply(1:K,
                        function(i)
                          trial_simul2(
                            N = Ns[i],
                            trt_names = network_settings$design[[i]],
                            delta = network_settings$delta[[i]],
                            deltasub = network_settings$subdelta[[i]],
                            mod_dist = network_settings$mod_prev[[i]],
                            sigma0 = sigmas[i],
                            seed = rseed+i
                          )$data |> 
                          add_column(study = as.character(i),
                                     .before = 1)
  ) |> 
    bind_rows()
  
  
  # estimand: average treatment effect. It is assumed that there is no prior knowledge of the effect modification (otherwise it would be most likely controlled for)
  
  # estimate: calculate average treatment effect for each study by adjusting for known prognostic factor BUT not for effect modifier
  
  # TWO-STAGE NMA (Alternative. MULTINMA)
  raw_lm <- split(ipd_network, ipd_network$study) |>
    map_df(
      \(df) lm(y~trt_name + x, data = df) |> 
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
  
  
  nma <- netmeta(TE, seTE, treat1, treat2, 
                 studlab, data = netdata)
  
  netmeta::netgraph(nma)
  
  # inconsistency test
  netmeta:::forest.netsplit(
    netmeta::netsplit(nma, show = "all")
  )
  
  cat("+++++++****** RESULTS *******")

  cat("\n True indirect effect in IMT: \n")
  print(imt$indir_eff)
  
  message("IMT: all head-to-head comparisons (LSM)")

  imtcntr <- emmeans::contrast(
    emmeans::emmeans(imtmod, ~trt_name) ,
    method = "pairwise", 
  )
  imtcontr <- imtcntr |> 
    as.data.frame() |> 
    dplyr::bind_cols(
      emmeans:::confint.emmGrid(imtcntr) |> 
        dplyr::select(dplyr::ends_with("CL"))
    ) |> 
    dplyr::rename(stat = t.ratio)
  
  print(
    imtcontr |> 
      dplyr::mutate(
        dplyr::across(
          dplyr::where(is.numeric), 
          \(x) round(x, 2)
        )
      )
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
    upper.CL = nma[[paste0("upper.", nmatype)]][upper.tri(nma[[paste0("TE.", nmatype)]])]
  )
  
  cat(
    " Number of studies: k =", nma$k, "\n",
    "Number of pairwise comparisons: m =", nma$m, "\n",
    "Number of treatments: n =", nma$n, "\n",
    "Number of designs: d =", nma$d, "\n"
  )
  
  print(nmacontr |> 
          dplyr::mutate(
            dplyr::across(
              dplyr::where(is.numeric), 
              \(x) round(x, 2)
            )
          )
  ) 
  

  return(
    list(
      IMT = imtmod,
      NMA = nma,
      data = list(imt = imt$data, 
                  ind_eff = imt$indir_eff,
                  ipd_net = ipd_network),
      est = imtcontr |> 
        tibble::add_column(
          evidence = "IMT"
        ) |> 
        dplyr::bind_rows(
          nmacontr |> 
            tibble::add_column(
              evidence = "NMA"
            )
        )
    )
  )
  
}
