# V-shaped network

source("./R/dependencies.R")
source("./R/trial_simul2.R")
source("./R/network_simul.R")


# GENERAL SETTINGS
##### IMT ####
# settings:
# N = 10000 large enough ideal trial
# average TRT effects = BA: -10; CA = -5 
# modifier = two balanced strata P(V1)=0.5, P(V2)=0.5, 
# treatment effect in V1 = 0
# treatment effect in V2 = BA: -20; CA: -10
# indirect effect BC, by consistency = BA-CA= -10 - (-5) = -5
## simulate network ######
## settings:
## N studies = 5
## designs = BA, CA (head-to-head only, no multiarm)
## each study has same prognostic and modifier distribution as IMT
## common effect = all effects are assumed fixed across different study populations
## sigma = residual variance can vary from study to study
# NOTE: we assume prognostic variable, X, is always adjusted for which already accounts for different X distributions across studies (10.1136/bmjebm-2022-111931)

ssizes <- list(small = 1,
               medium = 10,
               large = 100)

inconsistency <- list(
  none = rep(0.5, 5),
  mild = c(0.5, 0.8, 0.5, 0.5, 0.5),
  high = c(0.5, 0.8, 0.5, 0.8, 0.8) # more inconsistency in AC minority design
)



######################### raw results ####################

rawres1 <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      
      network_simul(
        network_settings = list(
          N = c(10, 50)*ssizes[[j]], 
          design = list(
            c("A", "B"), c("A", "B"), c("A", "B"),
            c("A", "C"), c("A", "C")
          ),
          delta = as.list(c(rep(-10, 5-2), rep(-5, 5-3))),
          subdelta = as.list(rep(0, 5)),
          mod_prev = as.list(inconsistency[[i]]),
          sigma = c(0.5, 3)
        ),
        nma_common = TRUE
      )
  )
); names(rawres1) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawres1[[i]]) <- names(ssizes) 


# EXTRACT CONTRASTS

res1 <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      
      rawres1[[i]][[j]]$est |> 
      tibble::add_column(
        inconsistency = i,
        samplesize = j
      )
  )
) |> 
  dplyr::bind_rows() |> 
  dplyr::mutate(
    samplesize = factor(samplesize, 
                        levels = c("small", "medium", "large")),
    inconsistency = factor(inconsistency,
                           levels = c("none", "mild", "high"))
  )


################

# plot 

res1 |> 
  dplyr::filter(evidence == "NMA") |> 
  ggplot(
    aes(x = samplesize, y = estimate, 
        colour = contrast, group = contrast) # , shape = evidence
  ) +
  # geom_pointrange(
  #   aes(ymin = lower.CL, ymax = upper.CL) ) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(vars(inconsistency)) +
  geom_hline(yintercept = 10, colour  ="red3", linetype = 2, alpha = 0.5) + 
  geom_hline(yintercept = 5, colour  ="red3", linetype = 2, alpha = 0.5)  +
  geom_hline(yintercept = -5, colour  ="red3", linetype = 2, alpha = 0.5) +
  geom_point(
    data = res1 |> 
      dplyr::filter(evidence == "IMT" & inconsistency == "none"),
    aes(x = samplesize, y = estimate), shape = 5, size = 4
  ) +
  ggplot2::labs(
    x = "Sample size", 
    title = "Estimation in a V-shaped treatment network for increasing inconsistency levels (Nr of studies: 5)", 
    caption = paste0(
      "Diamonds: IMT for large N; Dots: NMA (", 
      unique(res1$model),")"
    )
    ) 

#########################################################################################################
#########################################################################################################
# INTERPRETATION: 
# 1) inconsistency none. In medium sample size, bias in AB contrast is entirely due to random error (random precision and/or residual error). This alone is sufficient to cause inconsistency in BC contrast. As precision becomes uniformly high (large sample size), NMA estimates are unbiased.
# 2) inconsistent. Weight of direct/indirect evidence, which design is affected by bias and how many studies, all determine how bias propagates across the network.Focusing on large sample size, inconsistency affects only one study in majority design AB (mild scenario) also causing bias in BC contrast. Interestingly, by increasing bias for both AB and AC designs seems to cancel out bias in BC design. However, this should be regarded as a change outcome and it stresses how unpredictable bias propagation can be as bias and network size increases.  
#########################################################################################################
#########################################################################################################

###################################### 
# ####### BALANCE populations ########
# ---> assuming IPD is available
## approach 1: in a two-stage NMA, adjust for other known factor, V (only ATE).
## approach 2: re-balance studies using IPW methods (ATE and ATT)
## approach 3: run multi-level MNA 

# ---> assuming IPD is unavailable
## approach 1: in a two-stage NMA, adjust for other known factor, V.
## approach 2: run multi-level MNA using summary data (except one study; only ATT)
## approach 3: GCIPDR and IPW pseudodata

# EXTRACT NETWORK DATA

res1dat <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      
      rawres1[[i]][[j]]$data$ipd_net |> 
      tibble::add_column(
        inconsistency = i,
        samplesize = j
      )
  )
) |> 
  dplyr::bind_rows() |> 
  dplyr::mutate(
    samplesize = factor(samplesize, 
                        levels = c("small", "medium", "large")),
    inconsistency = factor(inconsistency,
                           levels = c("none", "mild", "high"))
  )

### EXCURSUS: 

# Marginal effect of TRT modifier, V. Average of: V2-V1 contrast in arm B (-20), V2-V1 contrast in arm C (-10), V2-V1 contrast in control arm A (0) -> (-20 - 10 + 0)/3 = -30/3 = -10 
summary(
  lm(y~trt_name + x + V, 
     data = rawres1$mild$large$data$imt)
)
# TRT-V interaction returns V2-V1 contrast in the respective arm, B and C (e.g., -20 and -10), TRT_B and TRT_C is respective effect in V1 (0), V2 is effect in control A (0)
summary(
  lm(y~trt_name*V + x, 
     data = rawres1$mild$large$data$imt)
)


# ATE estimand

prova <- res1dat |> 
  filter(
    inconsistency == "high",
    samplesize == "small"
  ) |> 
  run_two_stage_nma(
    study_level_model_formula = formula(y~trt_name + x + V)
  )


# ATT estimand

  