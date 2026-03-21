#' Simulate treatments network and corresponding ideal multiarm trial (IMT)
#' @param N - integer - size of IMT
#' @param delta - vector - TRT effects versus common reference in IMT. TRT names are automatically labelled with letters, A, B, C ...
#' @param mod_dist - vector - prevalence in TRT modifier strata in IMT
#' @param deltasub - vector - equal length of mod_dist; TRT effect in modifier's stratum in IMT  
#' @param network_settings - list - with elements: 
#' - K (number of studies) 
#' - N (vector - sample size in respective study)
#' - design (list - TRTs in each study (head-to-head or multiarm); first name is control arm)
#' - delta (list - TRT effect(s) in each study)
#' - subdelta (list - TRT effect(s) in study's modifier stratum)
#' - mode_prev (list - prevalences of modifier's strata in each study; must have equal lenght of corresponding list element in subdelta)
#' - sigma (vector - residual error in each study)
#' @description
#' Number of unique TRTs in network must match number of treatments in IMT, anchored to same control. Design must be any combination of the above treatment (head-to-head or multi-arm). Design in each study determines the network structure (V-shape, triangle, star, etc ...)
#' 


network_simul <- function(
    N = 10000,
    delta = c(-10, -5),
    mod_dist = 0.5,
    deltasub = 0,
    # more arguments from trial_simul2 possible (and in network setting)
    network_settings = list(
      K = 5,
      N = round(runif(K, min = 100, max = 500)), # by increasing N estim become more precise
      design = list(
        c("A", "B"), c("A", "B"), c("A", "B"),
        c("A", "C"), c("A", "C")
      ),
      delta = as.list(c(rep(-10, K-2), rep(-5, K-3))),
      subdelta = as.list(rep(0, K)),
      mod_prev = as.list(rep(0.5, K)),
      sigma = runif(K, 0.5, 3), # study specific precision
      seed = 12
    )
){
  
  
  #### IMT ####

  imt <- trial_simul2(
    N = N,
    delta = delta,
    mod_dist = mod_dist,
    deltasub = deltasub
  )
  
  # check
  
  imtmod <- lm(y~trt_name + x, data = imt$data)
  
  summary(imtmod)
  # check consistency
  update(imtmod, contrasts = list(
    trt_name = contr.treatment(
      n=LETTERS[1:3], base = 3)))
  
  
  ## simulate network
  
  # stack network data
  
  ipd_network <- lapply(1:5,
                        function(i)
                          trial_simul2(
                            N = settings$N[i],
                            trt_names = settings$design[[i]],
                            delta = settings$delta[[i]],
                            deltasub = settings$subdelta[[i]],
                            mod_dist = settings$mod_prev[[i]],
                            sigma0 = settings$sigma[i],
                            seed = seed+i
                          )$data |> 
                          add_column(study = as.character(i),
                                     .before = 1)
  ) |> 
    bind_rows()
  
  
  # estimand: average treatment effect. It is assumed that there is no prior knowledge of the effect modification (otherwise it would be most likely controlled for)
  
  # estimate: calculate average treatment effect for each study by adjusting for known prognostic factor BUT not for effect modifier
  
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
  netmeta::netleague(nma) 
  
  
  # inconsistency cannot be tested in open loop graph
  netmeta:::forest.netsplit(
    netmeta::netsplit(nma, show = "all")
  )
  
  print(nma)
  

  return(
    list(
      IMT = imtmod,
      NMA = nma
    )
  )
  
}
