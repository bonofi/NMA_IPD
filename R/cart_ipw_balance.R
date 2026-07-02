#' Performs IPW on CART data
#' 
#' 



cart_ipw_balance <- function(
    ipd_network,
    boot_iter = 100,
    seed = 7385
    
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
  
  
  
  
}



cart_ipw_balance(
  res1dat |> 
      filter(samplesize == "medium" & inconsistency == "high")
)