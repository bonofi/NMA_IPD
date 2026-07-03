#' Performs IPW on CART data
#' 
#' 



cart_ipw_balance <- function(
    ipd_network,
    modelformula = as.formula(study ~ x + V1 + V2),
    datalevel = c("ipd-agd", "agd"), # type of data available: mixed ipd-agd, only agd ...
    estimand = c("ATT", "ATE"),
    ref_study = "1",   # for ATT estimation
    stop_rule = "ks.mean",   # can be a vector
    n_trees = 5000,
    boot_iter = 100,
    seed = 7385,
    cores = 5
    
){
  
  
  browser()
  
 
  
  tictoc::tic()
  # synthetic data using original IPD 
  
  raw <- ipd_network |> 
    dplyr::select(
      y, study, V, trt_name, x) |> 
    synthpop::syn.strata(
      strata = "study",
      method = "parametric",  # alternative "cart"
      m=boot_iter, 
      seed = seed,
      minstratumsize = 10,
      ver
    )
  
  
  tictoc::toc()

  synthpop::summary.synds(raw)
  
  
  ## run IPW on synthetic data

  mirai::daemons(cores)
  
  # run IPW on pseudodata --> raw result
  tictoc::tic()
  rawipw <- lapply(raw$syn,
                   function(x)
                   {
                     out <- model.matrix(~ V - 1, data = x) |> 
                       as.data.frame() |> 
                       dplyr::bind_cols(
                         x
                       )
                       
                     
                   }
  ) |> 
    purrr::map(
      purrr::in_parallel(
        
        \(df) {
          
          ipwreg <- try(
            ipw_balance(
              df |> 
                # if datalevel = agd, it will be NULL
                dplyr::bind_rows(
                  refstudy
                ) |> 
                dplyr::arrange(study),
              model_formula = modelformula,
              estimand = estimand,
              stop_rule = stop_rule,
              ref_study = ref_study,
              n_trees = n_trees,
              print_diagnostics = FALSE
            ),
            silent = TRUE
          )
          
          if (class(ipwreg) == "try-error")
            return(
              list(
                est = data.frame(failed = ipwreg[1])
              )
            )
          else
            ipwreg
        }
        # pass environment objects to in_parallel
        , "%>%" = `%>%`,
        refstudy = refstudy,
        ipw_balance = ipw_balance,
        modelformula = modelformula,
        estimand = estimand,
        stop_rule = stop_rule,
        ref_study = ref_study,
        n_trees = n_trees
        ### end passing objects to in_parallel
      ) # end in_parallel
    )
  
  tictoc::toc()
  
  mirai::daemons(0)
  
  gc()
  
  
}



cart_ipw_balance(
  res1dat |> 
      filter(samplesize == "medium" & inconsistency == "high")
)