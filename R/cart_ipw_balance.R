#' Performs IPW on CART data
#' 
#' 


cart_ipw_balance <- function(
    ipd_network,
    modelformula = as.formula(study ~ x + V1 + V2),
    datalevel = c("ipd-agd", "agd"), # type of data available: mixed ipd-agd, only agd ...
    method = c("normrank", "", "cart", "sample", "cubertnorm" ),
    estimand = c("ATT", "ATE"),
    ref_study = "1",   # for ATT estimation
    stop_rule = "ks.mean",   # can be a vector
    n_trees = 5000,
    boot_iter = 100,
    seed = 7385,
    cores = 5,
    save_raw = FALSE
    
){
  
  
  datalevel <- match.arg(datalevel)
  estimand <- match.arg(estimand)
  
  studyid <- pickst <- unique(ipd_network$study)
  refstudy <- NULL
  
  # avoid simulation of reference study if pseudodata generation is internal
  if (datalevel == "ipd-agd")
  {
    pickst <- setdiff(studyid, ref_study)
    refstudy <- ipd_network |> 
      dplyr::filter(study == ref_study)
  }
 
  browser()
  # synthetic data using original IPD 
  
  raw <- ipd_network |> 
    dplyr::filter(study %in% pickst) |> 
    dplyr::select(
      y, study, V, trt_name, x) |> 
    synthpop::syn.strata(
      strata = "study",
      visit.sequence = c(3, 5, 4, 1, 2),
      method = method, #c("cart", "", "cart", "sample", "cart" ), #c("normrank", "", "cart", "sample", "cubertnorm" ),  
      m=boot_iter, 
      seed = seed,
      minstratumsize = 10
    )
  
  # todo: must iclude subjid!!
  
 # synthpop::summary.synds(raw)
 # synthpop::compare(
 #   raw,
 #   ipd_network |> 
 #     dplyr::filter(study %in% pickst) |> 
 #     dplyr::select(
 #       y, study, V, trt_name, x)
 #   )
  
  ## run IPW on synthetic data

  mirai::daemons(cores)
  
  # run IPW on pseudodata --> raw result
  tictoc::tic()
  
  rawipw <- lapply(raw$syn,
                   function(x)
                   {
                     out <- model.matrix(~V-1, data = x) |> 
                       as.data.frame() |> 
                       dplyr::bind_cols(
                         x |> 
                           dplyr::group_by(study) |> 
                           dplyr::mutate(
                             rownum = row_number(),
                             # need usubjid for compatibility with other utilities
                             usubjid = paste0(
                               study, "-", 
                               rownum),
                             .before = 1
                           ) |> 
                           dplyr::ungroup() |> 
                           dplyr::select(-rownum)
                       )
                     
                     colnames(out) <- stringr::str_replace_all(
                       colnames(out),
                       "Vlevel_", "V")
                     out
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
  
  
  # extract estimates and stack by boot iteration
  cleanipw <- 1:length(rawipw) |> 
    lapply(
      function(i)
        rawipw[[i]]$est |> 
        tibble::add_column(bootIter = i, 
                           .before = 1)
      
    ) |> 
    dplyr::bind_rows()
  
  # extract extra info
  extra <- data.frame(
    model = NA,
    evidence = "CART-IPW",
    estimand = estimand,
    level = toupper(datalevel),
    evidence2 = "Balanced"
  ) 
  
  # calculate summary over boot iteration
  
  summipw <- cleanipw |> 
    dplyr::select(-bootIter) |> 
    dplyr::group_by(
      contrast
    ) |> 
    dplyr::summarise(
      dplyr::across(
        # here names to be kept, e.g., estimate, lower, upper, etc ..
        .cols = dplyr::where(is.numeric),
        .fns = list(
          Mean = ~ mean(.x, na.rm = TRUE),
          SE = ~ sd(.x, na.rm = TRUE)
          # more ?
        ),
        .names = "{.col}_{.fn}"  # glue-style template
      )
    ) |> 
    dplyr::ungroup()
  
  meanest <- summipw |> 
    dplyr::select(
      contrast,
      dplyr::ends_with("_Mean")
    )
  pickm <- grepl("_Mean", 
                 names(meanest))
  names(meanest)[pickm] <- stringr::str_replace_all(
    names(meanest)[pickm], "_Mean", "")
  
  estse <- summipw |> 
    dplyr::select(
      contrast,
      dplyr::ends_with("_SE")
    )
  picks <- grepl("_SE", 
                 names(estse))
  names(estse)[picks] <- stringr::str_replace_all(
    names(estse)[picks], "_SE", "")
  
  
  return(
    list(
      corr_diagnostics_pseudodat = NULL,
      rawest = {
        if (save_raw)
          cleanipw
        else 
          NULL
      },
      est_se = estse |> 
        cbind(extra),
      est = meanest |> 
        cbind(extra)
    )
  )
  
  
  
}



test <- cart_ipw_balance(
  res1dat |> 
      filter(samplesize == "medium" & inconsistency == "high")
)