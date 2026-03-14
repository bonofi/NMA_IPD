#' Simulate Trial
#' @description
#' Simulate multi arm (arms >=2) trial from IPD generating mechanism
#' @description
#' Simulate trial with continuous primary outcome, Y, and residual unit-variance Normal error, e.g. cholesterol levels, in presence of prognostic variable, X, and treatment modifier strata, V.
#' @param N - scalar - overall sample size. It will be split across arms based on allocation ratio and number of treatment effect modifier strata.
#' @param mu0 - scalar - baseline average outcome Y.
#' @param beta - scalar - effect of prognostic variable X on outcome Y
#' @param delta - vector - of average treatment effects, one for each arm in the study. If delta is scalar, the study is head-to-head parallel.
#' @param deltasub - vector of treatment effects of K-1 ordered subgroups (effect modifications). The K-th effect is chosen such to average the prespecified delta effect.  
#' If delta has length > 1 (multiarm study) the deltasub vector contains the K-1 subgroup effect that are equal for all arms. For having fully customized arm-specific subgroup effects, let deltasub be a matrix with K-1 rows (effect modification) and for each column-arm (this assumes a common effect modifier for simplicity, if mod_dist is vector).
#' @param sigma0 - scalar - homoschedastic residual error
#' @param pt - vector of allocation ratios; should sum to 1; first element has to be allocation to reference arm. If scalar, is the allocation in a head-to-head trial, e.g., 1:1 -> 0.5.
#' @param fx - function to simulate prognostic variable X, e.g., age with mean 60 --> Gamma(60).
#' @param mod_dist - vector of prevalence of K-1 subgroups (strata percentages of effect modifier, V). K-th proportion is known from the K-1 ones. Could allow for treatment-specific effect modifier if put into matrix (Not tested).
#' @param no_modifier - Boolean - If TRUE, no modifier is used = Treatment effect is homogeneous across the whole population.


trial_simul2 <- function(N, delta, mu0, beta, deltasub = c(0, 10), sigma0 = 1, pt = NULL,
                        fx = function(x) rgamma(x, 60), mod_dist = c(0.2, 0.4),
                        seed = 5602783, trt_names = LETTERS[1:(length(delta) + 1)],
                        no_modifier = FALSE){
  
  ######### sanity checks ################

  # 1) deltasub and mod_dist conformity
  if (
    ifelse(is.null(dim(deltasub)),
           length(deltasub),
           dim(deltasub)[1]) != length(mod_dist)
  ) 
    stop("length of deltasub is not conform with length of mod_dist")
  
  # add last modifier stratum proportion
  mod_dist <- c(mod_dist, 1-sum(mod_dist))
  
  # 2) Allow for common K-1 effect modifications if not specified otherwise
  if (length(delta) > 1)
  {
    if (is.null(dim(deltasub))) 
      # if study is multi-arm (length delta > 1) but deltasub is not matrix, then you must assume same K-1 effect modifications across arms ....
    {
      message("Your study is multiarm but you are assuming a set of common K-1 effect modifications across arms. You can have fully customized effect modifications for each arm as corresponding new columns (arms) in a deltasub matrix.")
      # Assume K-1  common effect modifications if not stated otherwise 
      deltasub <- rep(deltasub, length(delta))
    } else if (dim(deltasub)[2] != length(delta)) {
      warning("Multiarm study: Number of columns of deltasub is not equal to number of arms !")
    }
    
  } 
  
  # 3) get allocation ratio if not given
  if (is.null(pt))
    #equal allocation
    pt <- rep(
      (N/(length(delta) + 1))/N,
      ifelse(
        length(delta) == 1,
        1,
        length(delta) + 1
      )
    )
  else if (length(pt) != (length(delta) + 1))
    warning("Allocation ratio not conform to number of arms!")
  
##### end sanity checks ################  
  
  set.seed(seed)
  
  res_err <- rnorm(N, sd = sigma0)
  
  # simulate prognostic variable
  x <- fx(N)
  
  # generate baseline population (target population) as the combination of average outcome mu0 and prognostic variable
  
  y0 <- mu0 + beta*x + res_err
  
  # distribution of effect modifier (not necessarily prognostic)
  
  z <- rmultinom(N, 1, prob = mod_dist) |>
    t()
  
  # split total trt effect into subgroup modified effects
  deltasub <- matrix(deltasub, ncol = length(delta))
  
  deltas <- rbind(
    deltasub,
    matrix( 
      (delta*dim(z)[2]) - apply(deltasub, 2, sum),
      nrow  =1)
  )
  
  # individual stratum-specific TRT effects before treatment assignment (at this stage any individual cannot be in all arms at the same time. TRT allocation is being performed below)   
  subdelta <- z%*%deltas
  

  # generate allocation list and 
  # Outcome expectation: linear relationship with effect modification
  if ( (length(delta) == length(pt)) == 1 )
  {
    # if head-to-head
    a <- rbinom(N, 1, pt)
    if (no_modifier)
      subdelta <- delta
  
    y <- y0 + subdelta*a
    trt_seq <- a
    aver_trt <- delta
    sub_delta <-  subdelta
    
  } else {
    
    # if multi-arm
    a <- rmultinom(N, 1, prob = pt) |>
      t()

    aver_trt <- a[, -1]%*%delta |> 
      as.vector()
    
    # drop reference (placebo) arm for calculating sub effects
    sub_delta <- (subdelta*a[, -1])%*%rep(1,length(delta)) |> 
      as.vector()
    
    if (no_modifier)
      sub_delta <- aver_trt
    y <- y0 + sub_delta 
    
    trt_seq <- (a%*%c(1:dim(a)[2]) - 1) |> 
      as.vector()
   
  }
     

  # simulate names
  
  trt_label <- data.frame(
    trt = sort(unique(trt_seq)),
    trt_name = sort(unique(trt_names))
  ) |> 
    dplyr::left_join(
      data.frame(
        trt = sort(unique(trt_seq))[-1],
        trt_eff = delta
      ),
      by = "trt"
    )
  

  # Get indirect effects
  indeff <- trt_label |> 
    na.omit() |> 
    dplyr::pull(trt_name) |> 
    combn(
      ifelse(length(delta) > 1, 
             2,
             1)
    ) |> 
    apply(2,
          function(cntr) {
            
            data.frame(
              contrast = paste(cntr, collapse = "-"),
              ind_eff = trt_label |> 
                dplyr::filter(
                  trt_name %in% cntr
                ) |> 
                dplyr::pull(trt_eff) |> 
                rev() |> 
                diff() |> 
                c(0) |> 
                sum()
            )
            
          },
          simplify = FALSE) |> 
    dplyr::bind_rows()
  
  
  out <- tibble::tibble(
    epsilon = res_err,
    mu0 = mu0,
    beta = beta,
    tot_delta = aver_trt,
    sub_delta = sub_delta,
    y = y,
    trt = trt_seq,
    x = x,
    V = paste(
      "level",
      z%*%c(1:dim(z)[2]),
      sep ="_"
    )
  ) |>
    cbind(
      as.data.frame(z)
    ) |> 
    dplyr::left_join(
      trt_label |> 
        dplyr::select(trt, trt_name),
      by = "trt"
    ) |> 
    tibble::add_column(
      ref_trt = paste0("seq:", 
                       trt_label$trt[1],
                       "; name:", 
                       trt_label$trt_name[1])
    )
  
  
  return(
    list(
      data = out,
      indir_eff = indeff
    )
  )
}

