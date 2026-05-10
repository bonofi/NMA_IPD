# V-shaped network

source("./R/dependencies.R")
source("./R/trial_simul2.R")
source("./R/network_simul.R")
source("./R/run_two_stage_nma.R")
source("./R/ipw_balance.R")
source("./R/multinma.R")
source("./R/gcipdr_ipw_balance.R")

# GENERAL SETTINGS
##### IMT ####
# settings:
# N = 10000 large enough ideal trial
# average TRT effects = BA: -10; CA = -5 
# modifier prevalence: we assume prevalence in general population, P(V1),
#  is equal to prevalence in trial = 0.5
# --> modifier = two balanced strata P(V1)=0.5, P(V2)=0.5, 
# treatment effect in V1 = 0
# treatment effect in V2 = BA: -20; CA: -10
# indirect effect BC, by consistency = BA-CA= -10 - (-5) = -5
# estimand: ATE = ATT
# reference trial = we assume reference trial is Nr 1 where P(V1) = 0.5
# ATE interpretation: if P(V1) not 0.5 in some trial, ATE = average effect across trials.
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
      unique(na.omit(res1$model)),")"
    )
    ) 

#########################################################################################################
#########################################################################################################
# INTERPRETATION: 
# 1) inconsistency none. In medium sample size, bias in AB contrast is entirely due to random error (random precision and/or residual error). This alone is sufficient to cause inconsistency in BC contrast. As precision becomes uniformly high (large sample size), NMA estimates are unbiased.
# 2) inconsistent. Weight of direct/indirect evidence, which design is affected by bias and how many studies, all determine how bias propagates across the network.Focusing on large sample size, inconsistency affects only one study in majority design AB (mild scenario) also causing bias in BC contrast. Interestingly, by increasing bias for both AB and AC designs seems to cancel out bias in BC design. However, this should be regarded as a change outcome and it stresses how unpredictable bias propagation can be as bias and network size increases.  
#########################################################################################################
#########################################################################################################



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



###################################### 
# ####### BALANCE populations ########
# #######################################
# ---> assuming IPD is available
# #######################################

## approach 1: in a two-stage NMA, adjust for other known factor, V (only ATE).
## approach 2: re-balance studies using IPW methods (ATE and ATT)
## approach 3: run ML-NMR 

############################################################
# ATE estimand: average effect across trials
############################################################


# raw balanced data: IPD available
# TODO: only balance for mild inconsistency as proof of concept ...?
# system.time(
#   
# rawbal1 <- split(res1dat, res1dat$inconsistency) |> 
#   purrr::map(
#     \(df1) split(df1, df1$samplesize) |> 
#       purrr::map(
#         \(df2){
#           
#           print(
#             paste0(
#               "inconsistency ", 
#               unique(df2$inconsistency), 
#               "; sample size ", unique(df2$samplesize)
#             )
#           )
#           
#           list(
#             "2sNMA" = run_two_stage_nma(
#               ipd_network = df2,
#               study_level_model_formula = formula(y~trt_name + x + V)
#             ),
#             "IPW" = ipw_balance(
#               ipd_network = df2,
#               model_formula = as.formula(study ~ x + V1 + V2),
#               estimand = "ATE",
#               stop_rule = "es.mean"
#             ),
#             "ML-NMR" = multinma(
#               ipd_network = df2,
#               modelformula = as.formula(~x + V),
#               datalevel = "ipd"
#             )
#           )
#         }
#           
#           
#       )
#   ) 
# )
# 
# names(rawbal1) <- names(inconsistency)
# for (i in names(inconsistency))
#   names(rawbal1[[i]]) <- names(ssizes) 

###### PARALLELIZE TASK 

simplan <- purrr::map(
  unique(res1dat$inconsistency),
  \(inc) purrr::map(
    unique(res1dat$samplesize),
    \(n) data.frame(
      inc = inc,
      n = n
    ) |> 
      dplyr::mutate(
        dplyr::across(
          everything(),
          ~as.character(.x)
        )
      )
  )
) |> 
  dplyr::bind_rows()
  

cores <- detectCores() - 1

future::plan(multisession, workers = cores)

tictoc::tic()

rawbal1 <- 1:dim(simplan)[1] |> 
  furrr::future_map(

  #  purrr::in_parallel(
      
      \(i){
        # set level
        inc <- simplan[i, "inc"]
        n <- simplan[i, "n"]   
        print(
          paste0(
            "inconsistency ", inc, 
            "; sample size ", n )
        )
        
        df2 <- res1dat |> 
          dplyr::filter(inconsistency == inc & 
                          samplesize == n)
        
        list(
          "2sNMA" = run_two_stage_nma(
            ipd_network = df2,
            study_level_model_formula = formula(y~trt_name + x + V),
            print_network = FALSE
          ),
          "IPW" = ipw_balance(
            ipd_network = df2,
            model_formula = as.formula(study ~ x + V1 + V2),
            estimand = "ATE",
            stop_rule = "es.mean",
            print_diagnostics = FALSE,
            save_raw = FALSE
          ),
          "ML-NMR" = multinma(
            ipd_network = df2,
            modelformula = as.formula(~ x + V),
            datalevel = "ipd",
            print_network = FALSE,
            save_raw = FALSE
          )
        )
      },
      .options = furrr_options(
          seed = TRUE,
          globals = c("estimand", "simplan", "res1dat",
                      "run_two_stage_nma", "ipw_balance", 
                      "multinma")) 
        
        #  )
        
      )

tictoc::toc()

names(rawbal1) <- apply(simplan, 1, \(x) paste(x, collapse = "-"))
 
future::plan(sequential)


# rawbal1 <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1.rds")

### extract contrasts after balancing

res1bal <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      # use this code
      # lapply(
      #   c("2sNMA", "IPW", "ML-NMR"),
      #   function(x)
      #     rawbal1[[i]][[j]][[x]]$est
      # ) |> 
      # dplyr::bind_rows()
      rawbal1[[i]][[j]]$`2sNMA`$table |>
      bind_rows(
        rawbal1[[i]][[j]]$IPW$est
      ) |> 
      bind_rows(
        rawbal1[[i]][[j]]$`ML-NMR`$est
      ) |> 
      tibble::as_tibble() |> 
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

# res1bal$evidence2[res1bal$evidence == "2-stage adjusted NMA"] <- "Balanced"

allres1 <- res1 |> 
  dplyr::bind_rows(
    res1bal
  )

#### plot


allres1 |> 
  dplyr::rename(Method = evidence) |> 
  dplyr::filter(evidence2 == "Balanced") |> 
  ggplot(
    aes(x = samplesize, y = estimate, 
        colour = contrast, group = contrast, shape = Method) 
  ) +
  geom_line(
    data = allres1 |> 
      dplyr::rename(Method = evidence) |>
      dplyr::filter(Method == "IPW")
  ) +
  geom_point(size = 2) +
  facet_wrap(vars(inconsistency)) +
  geom_hline(yintercept = 10, colour  ="red3", linetype = 2, alpha = 0.5) + 
  geom_hline(yintercept = 5, colour  ="green3", linetype = 2, alpha = 0.5)  +
  geom_hline(yintercept = -5, colour  ="skyblue4", linetype = 2, alpha = 0.5) +
  geom_point(
    data = res1 |> 
      dplyr::filter(evidence == "IMT" & inconsistency == "none"),
    aes(x = samplesize, y = estimate), shape = 5, size = 4
  ) +
  ggplot2::labs(
    x = "Sample size", 
    title = "ATE estimation in a V-shaped treatment network for increasing inconsistency levels (Nr of studies: 5) after balancing: IPD available", 
    caption = paste0(
      "Diamonds: IMT for large N; Dots: NMA (", 
      unique(na.omit(res1bbal$model)),")"
    )
  ) 

# INTERPRETATION: ATE is an average across trials that do not necessarily coincide with true population average iff trial average do not coincide with population average

#############################################################
# ATT estimand: effect in reference trial Nr 1
#############################################################
  

system.time(
  
  rawbal1b <- split(res1dat, res1dat$inconsistency) |> 
    purrr::map(
      \(df1) split(df1, df1$samplesize) |> 
        purrr::map(
          \(df2){
            
            print(
              paste0(
                "inconsistency ", 
                unique(df2$inconsistency), 
                "; sample size ", unique(df2$samplesize)
              )
            )
            
            list(
              # "2sNMA" = run_two_stage_nma(
              #   ipd_network = df2,
              #   study_level_model_formula = formula(y~trt_name + x + V)
              # ),
              "IPW" = ipw_balance(
                ipd_network = df2,
                model_formula = as.formula(study ~ x + V1 + V2),
                estimand = "ATT",
                stop_rule = "es.mean"
              )
              # ,
              # "ML-NMR" = multinma(
              #   ipd_network = df2,
              #   modelformula = as.formula(~x + V),
              #   datalevel = "ipd"
              # )
            )
          }
          
          
        )
    ) 
)

names(rawbal1b) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1b[[i]]) <- names(ssizes) 


##
res1bbal <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      # use this code
      lapply(
        c(
         # "2sNMA", 
          "IPW" 
         # "ML-NMR"
        ),
        function(x)
          rawbal1b[[i]][[j]][[x]]$est
      ) |>
      dplyr::bind_rows()|>
      tibble::as_tibble() |> 
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


allres1b <- allres1 |> 
  dplyr::bind_rows(
    res1bbal
  )



allres1b |> 
  dplyr::rename(Method = evidence) |> 
  dplyr::filter(evidence2 == "Balanced" & estimand == "ATT") |> 
  ggplot(
    aes(x = samplesize, y = estimate, 
        colour = contrast, group = contrast, shape = Method) 
  ) +
  geom_line(
    data = allres1b |> 
      dplyr::rename(Method = evidence) |>
      dplyr::filter(Method == "IPW" & evidence2 == "Balanced" & estimand == "ATT")
  ) +
  geom_point(size = 2) +
  facet_wrap(vars(inconsistency)) +
  geom_hline(yintercept = 10, colour  ="red3", linetype = 2, alpha = 0.5) + 
  geom_hline(yintercept = 5, colour  ="green3", linetype = 2, alpha = 0.5)  +
  geom_hline(yintercept = -5, colour  ="skyblue4", linetype = 2, alpha = 0.5) +
  geom_point(
    data = res1 |> 
      dplyr::filter(evidence == "IMT" & inconsistency == "none"),
    aes(x = samplesize, y = estimate), shape = 5, size = 4
  ) +
  ggplot2::labs(
    x = "Sample size", 
    title = "ATT estimation in a V-shaped treatment network for increasing inconsistency levels (Nr of studies: 5) after balancing: IPD available", 
    caption = paste0(
      "Diamonds: IMT for large N; Dots: NMA (", 
      unique(na.omit(res1bbal$model)),")"
    )
  ) 

########################################################################
# ---> assuming IPD is NOT available for all studies except study Nr 1
# ######################################################################

## approach 1: in a two-stage NMA, adjust for other known factor, V. SAME as ATE
## approach 2: run multi-level MNA using summary data (except one study) -> ATE
## approach 3: GCIPDR and IPW pseudodata -> ATT

system.time(
  
rawbal1c <- split(res1dat, res1dat$inconsistency) |> 
  purrr::map(
    \(df1) split(df1, df1$samplesize) |> 
      purrr::map(
        \(df2){
          
          print(
            paste0(
              "inconsistency ", 
              unique(df2$inconsistency), 
              "; sample size ", unique(df2$samplesize)
            )
          )
          
          list(
            
            "GC-IPW" = gcipdr_ipw_balance(
              ipd_network = df2,
              modelformula = as.formula(study ~ x + V1 + V2),
              estimand = "ATT",
              stop_rule = "es.mean",
              boot_iter = 100
            )

          )
        }
        
        
      )
  ) 
)

names(rawbal1c) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1c[[i]]) <- names(ssizes) 


res1cbal <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      # use this code
      lapply(
        c(
          "GC-IPW"
        ),
        function(x)
          rawbal1c[[i]][[j]][[x]]$est
      ) |>
      dplyr::bind_rows()|>
      tibble::as_tibble() |> 
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


allres1c <- allres1b |> 
  dplyr::bind_rows(
    res1cbal
  )


#todo: check if outer level of furrr can speed up module. Check if SI_only is better 

# rawbal1c <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1c.rds")


########################################################################
# ---> assuming IPD is NOT available for all studies
# ######################################################################

system.time(
  
  rawbal1d <- split(res1dat, res1dat$inconsistency) |> 
    purrr::map(
      \(df1) split(df1, df1$samplesize) |> 
        purrr::map(
          \(df2){
            
            print(
              paste0(
                "inconsistency ", 
                unique(df2$inconsistency), 
                "; sample size ", unique(df2$samplesize)
              )
            )
            
            list(
              
              "GC-IPW" = gcipdr_ipw_balance(
                ipd_network = df2,
                modelformula = as.formula(study ~ x + V1 + V2),
                estimand = "ATT",
                datalevel = "agd",
                stop_rule = "es.mean",
                boot_iter = 100
              )
              
            )
          }
          
          
        )
    ) 
)

names(rawbal1d) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1d[[i]]) <- names(ssizes) 

# rawbal1d <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1d.rds")
