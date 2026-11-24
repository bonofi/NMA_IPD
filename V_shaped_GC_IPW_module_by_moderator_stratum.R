# run GC method for ATT in level IPD-AGD and AGD by stratifying synthesis for moderator variable
# this module is aimed at showing that synthesizing data by subgroup can better capture the non-linear 
# relationship between the moderator and the outcome



## GC stratified by the mediator variable: proof of concept

########################################################
########################################################

future::plan(future::multisession, workers = detectCores() - 3)

B <- 300
system.time(
  
  rawGC3_strata <- lapply(
    c("11", "13"),
    \(m) split(res1dat, res1dat$inconsistency) |> 
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
              
              gc()
              # stratify by moderator
              out <- lapply(
                sort(unique(df2$V)), 
                function(l) {
                  
                  gc()
                  
                  do_gcipdr(
                    ipd_network = df2 |> 
                      dplyr::filter(V == l) |> 
                      dplyr::select(!starts_with("V")),
                    boot_iter = B,
                    method = m,
                    SI_k = c(50000, rep(10000, 4)), # only reference study 1 gets the optimal MC setting
                    only_SI = TRUE,
                    seed = 30697
                  )
                }
              )
              names(out) <- sort(unique(df2$V))
              
              # iterum componere strata in unum datum
              
              # matrix of moderator dummy variables
              Vdat <- diag(length(unique(df2$V)))
              rownames(Vdat) <- names(out)
              colnames(Vdat) <- stringr::str_replace_all(
                names(out), pattern = "level_", "V"
              )
              
              pseud <- lapply(
                1:B,
                function(b)
                  lapply(
                    names(out),
                    function(l)
                      out[[l]]$pseud[[b]] |> 
                      dplyr::select(
                        !dplyr::starts_with("V")
                      ) |> 
                      tibble::add_column(
                        !!!as.list(
                          Vdat[l, ]
                        )
                      )
                  ) |> 
                  dplyr::bind_rows() |> 
                  dplyr::arrange(study)
                
              ) 
              
              return(
                list(
                  raw = lapply(out, function(x) x$raw ),
                  pseud = pseud
                )
              )
            }     
          )
      ) 
    
  )
  
); names(rawGC3_strata) <- paste0("m_", c("11", "13")) 



set.seed(2607)
index <- sample(1:300, 100)

rawgc_att_strata <- lapply(
  names(rawGC3_strata),
  \(m)  lapply(
    c("ipd-agd", "agd"),
    \(dl) names(rawGC3_strata[[m]]) |> 
      purrr::map(
        \(i) names(rawGC3_strata[[m]][[i]]) |> 
          purrr::map(
            \(j){
              
              print(
                paste0(
                  "inconsistency ", i, 
                  "; sample size ", j
                )
              )
              gc()
              
              rlang::list2(
                
                "GC-IPW-{m}" := gcipdr_ipw_balance(
                  ipd_network = res1dat |> 
                    dplyr::filter(inconsistency == i & samplesize == j),
                  do_pseudodata = FALSE,
                  pseudodata = rawGC3_strata[[m]][[i]][[j]]$pseud[index],
                  modelformula = as.formula(study ~ x + V1 + V2),
                  estimand = "ATT",
                  datalevel = dl,
                  stop_rule = "es.mean"
                  
                )
                
              )
            }
            
          )
      ) 
  ) 
  
)
 
names(rawgc_att_strata) <- names(rawGC3_strata) 
  

for (m in names(rawgc_att_strata))
  names(rawgc_att_strata[[m]]) <- c("ipd-agd", "agd")

for (m in names(rawgc_att_strata))
  for (i in names(rawgc_att_strata[[m]])){
    for (j in names(inconsistency)){
      names(rawgc_att_strata[[m]][[i]]) <- names(inconsistency)
      names(rawgc_att_strata[[m]][[i]][[j]]) <- names(ssizes)
    }
  }

saveRDS(rawgc_att_strata, "./results/rawgc_att_strata.rds")

gc()



gc_att_strata <- lapply(
  names(rawgc_att_strata),
  \(m) lapply(
    names(rawgc_att_strata[[m]]),
    \(dl) lapply(
      names(inconsistency),
      \(i) lapply(
        names(ssizes), 
        \(j)
        
        rawgc_att_strata[[m]][[dl]][[i]][[j]][[paste0("GC-IPW-", m)]]$est |>
          tibble::as_tibble() |> 
          dplyr::mutate(
            inconsistency = i,
            samplesize = j,
            evidence = paste0("GC-IPW-strata-", m)
          )
      ) |> 
        dplyr::bind_rows()
    ) |> 
      dplyr::bind_rows() 
  ) |> 
    dplyr::bind_rows() 
) |> 
  dplyr::bind_rows() |> 
  dplyr::mutate(
    samplesize = factor(samplesize, 
                        levels = c("small", "medium", "large")),
    inconsistency = factor(inconsistency,
                           levels = c("none", "mild", "high"))
  )




#### some checks

rawgc_att_strata$m_11$`ipd-agd`$none$large$`GC-IPW-m_11`$rawest |> 
  ggplot(
    aes(
      x = estimate
    )
  ) + geom_histogram() + 
  facet_grid(
    rows = vars(contrast)
  )



rawgc_att_strata$m_13$`ipd-agd`$none$large$`GC-IPW-m_13`$rawest |> 
  ggplot(
    aes(
      x = estimate
    )
  ) + geom_histogram() + 
  facet_grid(
    rows = vars(contrast)
  )



# ############ OLD
# #### RUN GC-IPW IPD-AD (reference study nr 1 is in IPD format, all other are AD)
# 
# ## run IPW
# 
# system.time(
#   
#   rawbal1_attipdad_strata <- names(rawGC3_strata) |> 
#     purrr::map(
#       \(i) names(rawGC3_strata[[i]]) |> 
#         purrr::map(
#           \(j){
#             
#             print(
#               paste0(
#                 "inconsistency ", i, 
#                 "; sample size ", j
#               )
#             )
#             gc()
#             set.seed(2607)
#             list(
#               
#               "GC-IPW" = gcipdr_ipw_balance(
#                 ipd_network = rawGC3_strata[[i]][[j]]$pseud[sample(1:300, 100)],
#                 do_pseudodata = FALSE,
#                 modelformula = as.formula(study ~ x + V1 + V2),
#                 estimand = "ATT",
#                 stop_rule = "es.mean",
#                 cores = detectCores() - 2
#                 
#               )
#               
#             )
#           }
#           
#         )
#     ) 
# )
# 
# names(rawbal1_attipdad_strata) <- names(inconsistency)
# for (i in names(inconsistency))
#   names(rawbal1_attipdad_strata[[i]]) <- names(ssizes) 
# 
# gc()
# 
# rawbal1est_attipdad_strata <- lapply(
#   names(inconsistency),
#   function(i) lapply(
#     names(ssizes), 
#     function(j)
#       lapply(
#         c(
#           "GC-IPW-strata"
#         ),
#         function(x)
#           rawbal1_attipdad_strata[[i]][[j]][[x]]$est   
#       ) |>
#       dplyr::bind_rows()|>
#       tibble::as_tibble() |> 
#       tibble::add_column(
#         inconsistency = i,
#         samplesize = j
#       )
#   )
# ) |> 
#   dplyr::bind_rows() |> 
#   dplyr::mutate(
#     samplesize = factor(samplesize, 
#                         levels = c("small", "medium", "large")),
#     inconsistency = factor(inconsistency,
#                            levels = c("none", "mild", "high")),
#     evidence = "GC-IPW-strata"
#   ) 
# 
