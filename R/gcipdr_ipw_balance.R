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
      dplyr::select(where(is.numeric))
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
        NI_maxEval = 0)
    )
 
     )
  names(raw) <- unique(ipd_network$study)
  
    
  # pool pseudodata by study
  
  lapply(1:boot_iter, function(h)  # bootstrap's realization
    
    lapply(
      unique(ipd_network$study), 
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