#' perform IPW on GCIPD 
#' 
#' 
#' 


gcipd_ipw_balance(
  ipd_network,
  modelformula = as.formula(~x + V),   # ~(x + V)*.trt
  datalevel = c("ipd-agd", "agd"),
  estimand = c("ATT", "ATE")
  
){
  
  datalevel <- match.arg(datalevel)
  estimand <- match.arg(estimand)
  
  
}

