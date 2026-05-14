#' perform IPW on GCIPD 
#' 
#' 
#' 


gcipdr_ipw_balance <- function(
    ipd_network,
    modelformula = as.formula(study ~ x + V1 + V2),   
    datalevel = c("ipd-agd", "agd"), # type of data available: mixed ipd-agd, only agd ...
    estimand = c("ATT", "ATE"),
    ref_study = "1",   # for ATT estimation
    stop_rule = "ks.mean",   # can be a vector
    n_trees = 5000,
    method = c("3","4"),    # 3 = continuous Gamma, 4 = continuous Johnson
    boot_iter = 100,
    NI_maxEval = 10000,
    SI_k = 10000,
    only_SI = FALSE,
    seed = 49632,
    save_raw = TRUE,
    cores = detectCores() - 1 # ncores to use
    
){
  
  datalevel <- match.arg(datalevel)
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  
  
  studyid <- pickst <- unique(ipd_network$study)
  refstudy <- NULL
  
  if (datalevel == "ipd-agd")
  {
    pickst <- setdiff(studyid, ref_study)
    refstudy <- ipd_network |> 
      dplyr::filter(study == ref_study)
  }
  
  
  
  #set.seed(seed, "L'Ecuyer")
  pseudodata <- do_gcipdr(
    ipd_network |> 
      dplyr::filter(study %in% pickst),
    boot_iter = boot_iter,
    method = method,
    NI_maxEval = NI_maxEval,
    SI_k = SI_k,
    only_SI = only_SI,
    seed = seed,
    cores = cores
  )
  
  #future::plan(multisession, workers = cores)
  mirai::daemons(cores)
  
  # run IPW on pseudodata --> raw result
  tictoc::tic()
  rawipw <- pseudodata$pseud |> 
    #    furrr::future_map(      # Error: must change "future.globals.maxSize" options
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
  
  #future::plan(sequential)
  mirai::daemons(0)
  
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
  extra <- cleanipw |> 
    dplyr::distinct(model, evidence, estimand,
                    level, evidence2) |> 
    dplyr::mutate(
      evidence = "GC-IPW",
      level = toupper(datalevel)
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
  # eventually rearrange with pivot to have same format as ipw_balance$est ?
  
  meanest <- summipw |> 
    dplyr::select(
      contrast,
      dplyr::ends_with("_Mean")
    )
  names(meanest)[grepl("_Mean", 
                       names(meanest))] <- stringr::str_replace_all(
    names(meanest), "_Mean", ""
  )
  
  estse <- summipw |> 
    dplyr::select(
      contrast,
      dplyr::ends_with("_SE")
    )
  names(estse)[grepl("_SE", 
                     names(estse))] <- stringr::str_replace_all(
    names(estse), "_SE", ""
  )
  
  
  return(
    list(
      corr_diagnostics_pseudodat = pseudodata$raw |> 
        purrr::map(
          \(obj) obj$is.data.similar$lower.triangular.Rx
        ),
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


#' simulate pseudodata

do_gcipdr <- function(
    ipd_network,
    method = c("3","4"), # 3 = continuous Gamma, 4 = continuous Johnson
    boot_iter = 100,
    NI_maxEval = 200,
    SI_k = 10000,
    only_SI = FALSE,
    seed = 49632,
    cores = detectCores() - 1 # ncores to use
)
{
  
  method <- as.numeric(match.arg(method))
  
  trt_map <- ipd_network |> 
    dplyr::distinct(study, trt, trt_name)
  
  # prep data
  
  input <- lapply(
    split(ipd_network,
          unique(ipd_network$study)),
    function(x)
      x |> 
      dplyr::select(
        y, trt, x, starts_with("v")
      ) |> 
      dplyr::select(where(is.numeric)) |> 
      # drop redundant reference stratum level because 
      # it might cause trouble during optimization (corr values close to boundary)
      dplyr::select(!any_of("V1"))
  )
  
  # generate pseudodata. Output: list with boot repetition by study. Need to reorganize as list of pooled-by-study data repetitions  
  
  
  set.seed(seed, "L'Ecuyer") 
  
  future::plan(multisession, workers = cores)
  
  tictoc::tic()
  raw <- input |> 
    furrr::future_map(

        \(df) gcipdr::Simulate.many.datasets(
          list(df),
          H = boot_iter, 
          method = method, 
          checkdata = TRUE, 
          tabulate.similar.data = TRUE, 
          stochastic.integration = only_SI, # will override NI setup if TRUE
          NI_maxEval = NI_maxEval,
          SI_k = SI_k
        )[[1]],
        .options = furrr::furrr_options(
          seed = TRUE,
          globals = c("mclapply", "skewness", 
                      "adaptIntegrate", "rmvnorm"))
    )
  
  tictoc::toc()
  
  
  # check if any GC failed with numerical integration only
  fails <- which(
    unlist(
      lapply(raw, function(x) class(x))
    ) == "try-error"
  )
  
  # if some fails, rerun with MC integration
  if (length(fails) > 0){
    
    set.seed(seed, "L'Ecuyer") 
    
    tictoc::tic()
    raw[names(raw)[fails]] <- 
      names(raw)[fails] |> 
      furrr::future_map(
        
        \(i) gcipdr::Simulate.many.datasets(
          input[i],
          H = boot_iter, 
          method = method, 
          checkdata = TRUE, 
          tabulate.similar.data = TRUE, 
          stochastic.integration = TRUE, # will override NI setup if TRUE
          SI_k = SI_k
        )[[1]],
        .options = furrr::furrr_options(
          seed = TRUE,
          globals = c("mclapply", "skewness", 
                      "adaptIntegrate"))
      )
    
    tictoc::toc()
    
  }
  
  future::plan(sequential)
  
  # pool pseudodata by study
  
  out <- lapply(1:boot_iter, 
                function(h)  # bootstrap's realization
                  
                  lapply(
                    names(raw), 
                    function(j) ##  row-bind by study
                      
                      as.data.frame(
                        raw[[j]]$similar.data[[h]]
                      ) |> 
                      tibble::add_column(
                        study = j,
                        # need usubjid for compatibility with other utilities
                        usubjid = paste0(
                          j, "-", 
                          1:dim(raw[[j]]$similar.data[[h]])[1]),
                        .before = 1
                      ) |> 
                      # re-merge trt LABEL by study
                      dplyr::left_join(
                        trt_map |> 
                          dplyr::filter(
                            study == j
                          ),
                        by = c("study", "trt")
                      ) %>%
                      # must resort to colSums because rowwise is extremely slow !!!
                      dplyr::mutate(
                        # collect all V strata that are not the reference one
                        notV1 = rowSums(
                          . |> 
                            tibble::add_column(
                              # to be able to use rowSums in case of 1-dim V 
                              V0 = NA 
                            ) |> 
                            dplyr::select(
                              starts_with("V")),
                          na.rm = TRUE)) |>
                      dplyr::mutate(V1 = 1-notV1) |>
                      # drop notV1
                      dplyr::select(-notV1)
                    
                  ) |> 
                  dplyr::bind_rows()
                
  )
  
  
  return(
    list(
      raw = raw,
      pseud = out
    )
  )
  
  
}


# 
# prova <- gcipdr_ipw_balance(
#    res1dat |> 
#      filter(samplesize == "small" & inconsistency == "high"),
#    datalevel = "agd",
#    estimand = "ATT"
#    )