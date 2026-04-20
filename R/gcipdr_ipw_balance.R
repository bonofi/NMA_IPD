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
  
  # generate pseudodata

  set.seed(seed, "L'Ecuyer") # delete this line to assess stability
  
  print( 
    system.time(
      
      out <-  Simulate.many.datasets(
        split(ipd_network,
              unique(ipd_network$study)), 
        H = boot_iter, 
        method = 3,
        checkdata = TRUE, 
        tabulate.similar.data = TRUE, 
        NI_maxEval = 0)[[1]]
    )
 
     )
  )
    
    
  
  # pool pseudodata
  PoolArtifData <-function(istsimul.by.region){
    
    out <- lapply(1:2, function(i) lapply(1:100, function(h){
      
      merged <- do.call("rbind",
                        lapply(1:length(hubs), function(j){  ##  bind by country
                          
                          data <- as.data.frame(istsimul.by.region[[i]][[j]]$similar.data[[h]])  # matrix simul
                          data$COUNTRY <- hubs[j]
                          return(data) 
                        } )  )
      
      return(merged)
    })    ); names(out) <- paste0("meth_", 3:4)
    
    return(out)
  }
  
}




do_gcipdr(res1dat)