# run GC method for ATT in level IPD-AGD and AGD by stratifying synthesis for moderator variable
# this module is aimed at showing that synthesizing data by subgroup can better capture the non-linear 
# relationship between the moderator and the outcome



## GC stratified by the mediator variable: proof of concept

########################################################
########################################################


system.time(
  
  rawGC3_strata <- split(res1dat, res1dat$inconsistency) |> 
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
            
            B <- 300
            # stratify by moderator
            out <- lapply(
              c("level_1", "level_2"), # replace with levels
              function(l) {
                
                gc()
                
                do_gcipdr(
                  ipd_network = df2 |> 
                    dplyr::filter(V == l) |> 
                    dplyr::select(!starts_with("V")),
                  boot_iter = B,
                  method = "3",
                  SI_k = 30000,
                  only_SI = TRUE,
                  seed = 30697,
                  cores = detectCores()-2
                )
              }
            )
            names(out) <- c("level_1", "level_2")
            
            # conjunct strata into one single data set
            
            pseud <- lapply(
              1:B,
              function(b)
                out$level_1$pseud[[b]] |> 
                dplyr::select(
                  !dplyr::starts_with("V")
                ) |> 
                tibble::add_column(V1 = 1, V2 = 0) |>  # replace with bang-bang and diagonal 0-1 V matrix ?
                dplyr::bind_rows(
                  out$level_2$pseud[[b]] |> 
                    dplyr::select(
                      !dplyr::starts_with("V")
                    ) |> 
                    tibble::add_column(V1 = 0, V2 = 1)
                ) |> 
                arrange(study)
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





# ###############################################
# system.time(
#   
#   rawGC3_strata <- lapply(
#     c("level_1", "level_2"),
#     function(l)
#       
#       split(
#         res1dat |> 
#           dplyr::filter(V == l) |> 
#           dplyr::select(!starts_with("V")), 
#         res1dat$inconsistency[res1dat$V == l]
#       ) |> 
#       purrr::map(
#         \(df1) split(df1, df1$samplesize) |> 
#           purrr::map(
#             \(df2){
#               
#               print(
#                 paste0(
#                   "inconsistency ", 
#                   unique(df2$inconsistency), 
#                   "; sample size ", unique(df2$samplesize)
#                 )
#               )
#               
#               gc()
#               
#               do_gcipdr(
#                 ipd_network = df2,
#                 boot_iter = B,
#                 method = "3",
#                 SI_k = 30000,
#                 only_SI = TRUE,
#                 seed = 30697,
#                 cores = detectCores()-2
#               )
#               
#             }
#             
#             
#           )
#       ) 
#   )
#   
# ) 
# 
# names(rawGC3_strata) <- c("level_1", "level_2")
# 
# # regroup by mediator strata
# 
# rawGC3 <- names(rawGC3_strata[[1]]) |> 
#   purrr::map(
#       \(i) names(rawGC3_strata[[1]][[i]]) |> 
#         purrr::map(
#           \(j) 1:B |> 
#             purrr::map(
#               \(b) rawGC3_strata[["level_1"]][[i]][[j]][[b]] |> 
#                 dplyr::select(!starts_with("V")) |> 
#                 tibble::add_column(
#                   V1 = 1, V2 = 0
#                 ) |> 
#                 dplyr::bind_rows(
#                   rawGC3_strata[["level_2"]][[i]][[j]][[b]] |> 
#                     dplyr::select(!starts_with("V")) |> 
#                     tibble::add_column(
#                       V1 = 0, V2 = 1
#                     )
#                     
#                 )
#             )
#         )
#     )


#### RUN GC-IPW IPD-AD (reference study nr 1 is in IPD format, all other are AD)



## run IPW

system.time(
  
  rawbal1_attipdad <- names(rawGC) |> 
    purrr::map(
      \(i) names(rawGC[[i]]) |> 
        purrr::map(
          \(j){
            
            print(
              paste0(
                "inconsistency ", i, 
                "; sample size ", j
              )
            )
            gc()
            set.seed(2607)
            list(
              
              "GC-IPW" = gcipdr_ipw_balance(
                ipd_network = rawGC[[i]][[j]]$pseud[sample(1:300, 100)],
                do_pseudodata = FALSE,
                modelformula = as.formula(study ~ x + V1 + V2),
                estimand = "ATT",
                stop_rule = "es.mean",
                cores = detectCores() - 2
                
              )
              
            )
          }
          
        )
    ) 
)

names(rawbal1_attipdad) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1_attipdad[[i]]) <- names(ssizes) 

gc()

res1_attipdad <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      lapply(
        c(
          "GC-IPW"
        ),
        function(x)
          rawbal1_attipdad[[i]][[j]][[x]]$est   
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
                           levels = c("none", "mild", "high")),
    evidence = "GC-IPW"
  ) 

