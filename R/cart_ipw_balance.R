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
    
  
  tictoc::toc()
  
  
  prova <-  synthpop::syn(ipd_network |> 
                            
                            dplyr::filter(study == "1") |> 
                            select(y, x, V1, trt) ,
                          
                          m= 1)
  
  
}



cart_ipw_balance(
  res1dat |> 
      filter(samplesize == "medium" & inconsistency == "high")
)