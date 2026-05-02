#' perform IPW on GCIPD 
#' 
#' 
#' 


gcipdr_ipw_balance <- function(
  ipd_network,
  modelformula = as.formula(~x + V),   # ~(x + V)*.trt
  datalevel = c("ipd-agd", "agd"),
  estimand = c("ATT", "ATE")
  
){
  
  datalevel <- match.arg(datalevel)
  estimand <- match.arg(estimand)
  
  
}


#' simulate pseudodata

do_gcipdr <- function(
    ipd_network,
    boot_iter = 100,
    seed = 49632
    )
{
  
  browser()
  
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
        NI_maxEval = 0)
    )
 
     )
  names(raw) <- unique(ipd_network$study)
  
  # check if any GC failed with numericla integration only
  fails <- which(
    unlist(
      lapply(raw, function(x) class(x))
    ) == "try-error"
  )
    

  # if some fails, rerun with MC integration
  if (length(fails) > 0){
    for (i in fails)
      print( 
        system.time(
          raw[[i]] <- gcipdr::Simulate.many.datasets(
            input[i],
            H = boot_iter, 
            method = 3, # NORTA with Gamma marginals and Pearson corr
            checkdata = TRUE, 
            tabulate.similar.data = TRUE,
            stochastic.integration = TRUE,
            SI_k = 20000

            )
        )
        
      )
        
  }
      
  
    
    
  # pool pseudodata by study
  
  out <- lapply(1:boot_iter, function(h)  # bootstrap's realization
    
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
          
        )
    ) |> 
      dplyr::bind_rows()
    
    
  )
  
}
  
  





do_gcipdr(res1dat)