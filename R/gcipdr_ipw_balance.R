#' perform IPW on GCIPD 
#' 
#' 
#' 


gcipdr_ipw_balance <- function(
    ipd_network,
    do_pseudodata = TRUE, # if false additionally give a list of pseudodata 
    pseudodata = NULL,
    modelformula = as.formula(study ~ x + V1 + V2),   
    datalevel = c("ipd-agd", "agd"), # type of data available: mixed ipd-agd, only agd ...
    estimand = c("ATT", "ATE"),
    ref_study = "1",   # for ATT estimation
    stop_rule = "ks.mean",   # can be a vector
    n_trees = 5000,
    method = c("3","4","13"),    # 3 = continuous Gamma, 4 = continuous Johnson, 13 = user-defined marginals 
    boot_iter = 100,
    NI_maxEval = 10000, # or a named list "1" = 10000, "2" = 30000, ...
    SI_k = 10000, # or a named list "1" = 10000, "2" = 30000, ...
    only_SI = FALSE,
    seed = 49632,
    save_raw = TRUE,
    cores = detectCores() - 1, # ncores to use
    drop_ref_V = "V1",
    kruskal_use = FALSE,
    nr_future = FALSE
    
){
  
  datalevel <- match.arg(datalevel)
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  
  
  studyid <- pickst <- unique(ipd_network$study)
  refstudy <- NULL
  
  # avoid simulation of reference study if scheme is ipd-agd
  if (datalevel == "ipd-agd")
  {
    pickst <- setdiff(studyid, ref_study)
    refstudy <- ipd_network |> 
      dplyr::filter(study == ref_study)
  }
  
  if (do_pseudodata){
    pseudodata <- do_gcipdr(
      ipd_network |> 
        dplyr::filter(study %in% pickst),
      boot_iter = boot_iter,
      method = method,
      NI_maxEval = NI_maxEval,
      SI_k = SI_k,
      only_SI = only_SI,
      seed = seed,
      cores = cores,
      drop_ref_V = drop_ref_V,
      kruskal_use = kruskal_use,
      nr_future = nr_future,
    )
    
    gc()
    
  } else {
    if (is.list(pseudodata))
      pseudodata <- switch(
        datalevel,
        "agd" = list(pseud = pseudodata, raw = NULL),
        "ipd-agd" = list(
          # drop ref pseudostudy if ipd-agd
          pseud = lapply(
            pseudodata,
            \(x) x |> 
              dplyr::filter(study %in% pickst)
          ), 
          raw = NULL)  
      )
    else
      stop("do_pseudodata = FALSE but no pseudodata list provided")
  }

  future::plan(multisession, workers = cores)
  # mirai::daemons(cores)
  # 
  # #browser()
  # # run IPW on pseudodata --> raw result
  # tictoc::tic()
  # rawipw <- pseudodata$pseud |> 
  #   #    furrr::future_map(      # Error: must change "future.globals.maxSize" options
  #   purrr::map(
  #     purrr::in_parallel(
  #       
  #       \(df) {
  #         
  #         ipwreg <- try(
  #           ipw_balance(
  #             df |> 
  #               # if datalevel = agd, it will be NULL
  #               dplyr::bind_rows(
  #                 refstudy
  #               ) |> 
  #               dplyr::arrange(study),
  #             model_formula = modelformula,
  #             estimand = estimand,
  #             stop_rule = stop_rule,
  #             ref_study = ref_study,
  #             n_trees = n_trees,
  #             print_diagnostics = FALSE
  #           ),
  #           silent = TRUE
  #         )
  #         
  #         if (class(ipwreg) == "try-error")
  #           return(
  #             list(
  #               est = data.frame(failed = ipwreg[1])
  #             )
  #           )
  #         else
  #           ipwreg
  #       }
  #       # pass environment objects to in_parallel
  #       , "%>%" = `%>%`,
  #       refstudy = refstudy,
  #       ipw_balance = ipw_balance,
  #       modelformula = modelformula,
  #       estimand = estimand,
  #       stop_rule = stop_rule,
  #       ref_study = ref_study,
  #       n_trees = n_trees
  #       ### end passing objects to in_parallel
  #     ) # end in_parallel
  #   )
  # 
  # tictoc::toc()
  # 
  # #future::plan(sequential)
  # mirai::daemons(0)
  # 
  
  rawipw <- future.apply::future_lapply(
    pseudodata$pseud, 
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
    }, 
    future.seed = TRUE)
  
  
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
    evidence = "GC-IPW",
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
    method = c("3","4","11","13", "15"), # 3 = continuous Gamma, 4 = continuous Johnson, 13 = user-defined marginals
    boot_iter = 100,
    NI_maxEval = 200, # or a named list "1" = 100, "2" = 300 ...
    SI_k = 10000,  # or a named list "1" = 10000, "2" = 30000, ...
    only_SI = FALSE,
    seed = 49632,
    cores = detectCores() - 1, # ncores to use,
    drop_ref_V = "V1", # drop reference level of moderator variable -> keeping in increases instability due to perfect correlation
    kruskal_use = FALSE,
    nr_future = FALSE,
    interaction_only = FALSE # if TRUE, keeps only "inter" = trt*V and dropr trt and V

)
{
  
  method <- as.numeric(match.arg(method))
  
  trt_map <- ipd_network |> 
    dplyr::distinct(study, trt, trt_name)
  
  if (interaction_only)
    ipd_network <- ipd_network |> 
    dplyr::select(-trt, -dplyr::starts_with("V"))
  
  # prep data
  input <- lapply(
    split(ipd_network,
          ipd_network$study),
    function(x)
      x |> 
      dplyr::select(
        y, dplyr::starts_with("trt"), x, 
        dplyr::starts_with("v"), 
        dplyr::starts_with("inter")  # automatically loads any interaction term (by convention with prefix "inter")
      ) |> 
      dplyr::select(where(is.numeric)) |> 
      # drop redundant reference stratum level because 
      # it might cause trouble during optimization (corr values close to boundary)
      dplyr::select(!any_of(drop_ref_V)) 
  )

  # calculate marginal moments for custom marginals if method = "11", 13", "15"
  marginals <- NULL
  if (
    any(
      method %in% c("11", "13", "15")
    )
  )
    # internally user defined marginals
    marginals <- lapply(
      input,
      \(df) {
        moms <- df |> 
          dplyr::summarise(
            dplyr::across(
              .cols = dplyr::everything(),
              .fns = list(mean = mean, stdr = sd),
              .names = "{.col}_{.fn}",
              na.rm = TRUE
            )
          ) 
        # automatically detect all V moderators names from df
        Vnames <- grep("V", colnames(df), value = TRUE)
        out <- c() # empty list
        notV <- list(
          y = \(x) qnorm(x, mean = moms$y_mean, sd = moms$y_stdr),
          x = \(x) gcipdr:::qgmom(x, mean = moms$x_mean, sd = moms$x_stdr ),
          trt = \(x) qbinom(x, 1, prob = moms$trt_mean)
          # BEWARE: number of V moderators is hard-coded here but should automatically
          # adapt --> DONE
          #V2 = \(x) qbinom(x, 1, prob = moms$V2_mean)
        )
        
        allV <- lapply(
          Vnames,
          \(y) \(x) qbinom(x, 1, prob = moms[[paste0(y, "_mean")]]) # print(moms[[paste0(y, "_mean")]]) # test
        ); names(allV) <- Vnames
        
        # output and reorder list by df name
        c(out, notV, allV)[colnames(df)]
        
      }
      
    )
  
  
  # empty data container in case GC fails
  mockl <- rep(NA, dim(input[[1]])[2]) |> 
    as.list()
  names(mockl) <- colnames(input[[1]])
  mockd <- as.data.frame(mockl)
  
  
  # allow for study-specific gicpdr generative setting
  if (length(SI_k)  == 1)
  {
    SI_k <- rep(SI_k, length(input))
    names(SI_k) <- names(input)
  } else {
    if (
      !all(sort(names(SI_k)) == sort(names(input)))
    )
      stop("each element of the SI_k vector must is named after its corresponding study")
  }
  if (length(NI_maxEval)  == 1)
  {
    NI_maxEval <- rep(NI_maxEval, length(input))
    names(NI_maxEval) <- names(input)
  } else {
    if (
      !all(sort(names(NI_maxEval)) == sort(names(input)))
    )
      stop("each element of the NI_maxEval vector must is named after its corresponding study")
  }
  
    
  # generate pseudodata. Output: list with boot repetition by study. Need to reorganize as list of pooled-by-study data repetitions  
  
  set.seed(seed, "L'Ecuyer") 
  
  future::plan(multisession, workers = cores)
  
  tictoc::tic()
  raw <- names(input) |> 
    furrr::future_map(

        \(i) try(
          gcipdr::Simulate.data.given.IPD(
            input[[i]],
            H = boot_iter, 
            method = method, 
            checkdata = TRUE, 
            tabulate.similar.data = TRUE, 
            stochastic.integration = only_SI, # will override NI setup if TRUE
            NI_maxEval = NI_maxEval[i],
            SI_k = SI_k[i],
            user_defined_marginals = marginals[[i]],
            kruskal_use = kruskal_use,
            newtrap_parallel = nr_future
          ),
          silent = TRUE
        ),
        .options = furrr::furrr_options(
          seed = TRUE,
          globals = c("mclapply", "skewness", 
                      "adaptIntegrate", "rmvnorm",
                      "FitJohnsonDistribution"))
    )
  names(raw) <- names(input)
  
  tictoc::toc()
  
  
  # check if any GC failed with numerical integration only
  fails <- which(
    unlist(
      lapply(raw, function(x) class(x))
    ) == "try-error"
  )
  
  
  # if some fails, rerun with MC integration
  if (length(fails) > 0 & only_SI == FALSE){
    
    set.seed(seed, "L'Ecuyer") 
    
    tictoc::tic()
    raw[names(raw)[fails]] <- 
      names(raw)[fails] |> 
      furrr::future_map(
        
        \(i) try(
          gcipdr::Simulate.data.given.IPD(
            input[[i]],
            H = boot_iter, 
            method = method, 
            checkdata = TRUE, 
            tabulate.similar.data = TRUE, 
            stochastic.integration = TRUE, # will override NI setup if TRUE
            SI_k = SI_k[i],
            user_defined_marginals = marginals[[i]],
            kruskal_use = kruskal_use,
            newtrap_parallel = nr_future
          ),
          silent = TRUE
        ),
        .options = furrr::furrr_options(
          seed = TRUE,
          globals = c("mclapply", "skewness", 
                      "adaptIntegrate", "rmvnorm",
                      "FitJohnsonDistribution"))
      )
    
    tictoc::toc()
    
  }
  
  future::plan(sequential)
  
  
  # re-check if any GC failed with numerical integration only
  fails <- which(
    unlist(
      lapply(raw, function(x) class(x))
    ) == "try-error"
  )
  
  if (length(fails) > 0)
  {
    
    message(paste0("NO DATA PRODUCED FOR STUDY: ", 
                   paste(names(raw)[fails],
                         collapse = ", ")
    ))
    
    # Set mock shells for failed instances
    
    raw[names(raw)[fails]] <- names(raw)[fails] |>
      purrr::map(
        \(j) list(
          similar.data = lapply(
            1:boot_iter, function(b)
              mockd 
          )
        )
      )
    
  }
  
  
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
                        ) %>% 
                        {
                          if (interaction_only)
                            .
                          else
                            # re-merge trt LABEL by study
                            dplyr::left_join(
                              .,
                              trt_map |> 
                                dplyr::filter(
                                  study == j
                                ),
                              by = c("study", "trt")
                            )
                          
                        } %>%
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
                        # reintroduce complement V level
                        dplyr::mutate(
                          "{drop_ref_V}" := 1-notV1
                        ) |>
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