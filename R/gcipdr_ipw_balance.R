#' perform IPW on GCIPD 
#' 
#' 
#' 


gcipdr_ipw_balance(
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
    ipd_network
    )
{
  
  
  
  # generate pseudodata

  hubs <- levels(ist$COUNTRY)
  
  
  ######### LOAD COPULA DATA ##########
  
  
  seed <- 49632
  
  #
  istsimul.by.region <- lapply( 3:4, function(i)
    lapply(hubs, function(j){  
      
      data <- ist[ist$COUNTRY == j, -9]
      
      jiseed <- as.integer(paste(which(hubs == j), seed, i, sep = ""))
      set.seed( jiseed, "L'Ecuyer") # delete this line to assess stability
      print( system.time(
        artificial_data_object <-  Simulate.many.datasets(list(data), H = NULL, i,
                                                          checkdata = TRUE, tabulate.similar.data = TRUE, NI_maxEval = 0 )
      ))
      return( artificial_data_object[[1]])
    }  )         ); names(istsimul.by.region) <- paste0("meth_", 3:4)
  
  
  for(i in 1:2) names(istsimul.by.region[[i]]) <- hubs
  
  
  
  
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