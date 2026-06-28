#' Performs IPW on CART data
#' 
#' 



cart_ipw_balance <- function(
    ipd_network,
    boot_iter = 100,
    seed = 7385,
    cores = 5
    
){
  
  
  browser()
  
  future::plan(multisession, workers = cores)
  
  tictoc::tic()
  # synthetic data using original IPD 
  
  raw <- ipd_network |> 
    dplyr::select(
      y, study, V, trt_name, x) |> 
 synthpop::syn.strata(
        strata = "study",
        m=boot_iter, 
        seed = seed,
        minstratumsize = 10
        )
    )
  
  tictoc::toc()
  
  
  
}



cart_ipw_balance(
  res1dat |> 
      filter(samplesize == "small" & inconsistency == "high")
)