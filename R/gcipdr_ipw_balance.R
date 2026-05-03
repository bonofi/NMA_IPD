#' perform IPW on GCIPD 
#' 
#' 
#' 


gcipdr_ipw_balance <- function(
    ipd_network,
    modelformula = as.formula(~x + V1 + V2),   
    datalevel = c("ipd-agd", "agd"),
    estimand = c("ATT", "ATE"),
    boot_iter = 100,
    NI_maxEval = 200,
    SI_k = 10000,
    seed = 49632
    
){
  
  datalevel <- match.arg(datalevel)
  estimand <- match.arg(estimand)
  
  pseudodata <- do_gcipdr(
    ipd_network,
    boot_iter = boot_iter,
    NI_maxEval = NI_maxEval,
    SI_k = SI_k,
    seed = seed
  )
  
  
  
  
}


#' simulate pseudodata

do_gcipdr <- function(
    ipd_network,
    boot_iter = 100,
    NI_maxEval = 200,
    SI_k = 10000,
    seed = 49632
)
{
  
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
  
  print( 
    system.time(
      
      raw <-  gcipdr::Simulate.many.datasets(
        input,
        H = boot_iter, 
        method = 3, # NORTA with Gamma marginals and Pearson corr
        checkdata = TRUE, 
        tabulate.similar.data = TRUE, 
        stochastic.integration = FALSE,
        NI_maxEval = NI_maxEval)
    )
    
  )
  names(raw) <- unique(ipd_network$study)
  
  # check if any GC failed with numerical integration only
  fails <- which(
    unlist(
      lapply(raw, function(x) class(x))
    ) == "try-error"
  )
  
  
  # if some fails, rerun with MC integration
  if (length(fails) > 0){
    set.seed(seed, "L'Ecuyer") 
    
    for (i in fails)
      print( 
        system.time(
          raw[i] <- gcipdr::Simulate.many.datasets(
            input[i],
            H = boot_iter, 
            method = 3, # NORTA with Gamma marginals and Pearson corr
            checkdata = TRUE, 
            tabulate.similar.data = TRUE,
            stochastic.integration = TRUE,
            SI_k = SI_k
            
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
                          study = j
                        ) |> 
                        # re-merge trt LABEL by study
                        dplyr::left_join(
                          trt_map |> 
                            dplyr::filter(
                              study == j
                            ),
                          by = c("study", "trt")
                        ) %>%
                        dplyr::mutate(
                          # collect all V strata that are not the reference one
                          notV1 = rowSums(
                            . |> 
                              # to be able to use rowSums in case of 1-dim V
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
  
  
  return(out)
  
  
}



#do_gcipdr(res1dat)