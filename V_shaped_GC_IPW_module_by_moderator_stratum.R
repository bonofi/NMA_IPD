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
              sort(unique(df2$V)), 
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
            names(out) <- sort(unique(df2$V))
            
            # conjunct strata into one single data set
            
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

#### RUN GC-IPW IPD-AD (reference study nr 1 is in IPD format, all other are AD)

## run IPW

system.time(
  
  rawbal1_attipdad_strata <- names(rawGC3_strata) |> 
    purrr::map(
      \(i) names(rawGC3_strata[[i]]) |> 
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
              
              "GC-IPW-strata" = gcipdr_ipw_balance(
                ipd_network = rawGC3_strata[[i]][[j]]$pseud[sample(1:300, 100)],
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

names(rawbal1_attipdad_strata) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1_attipdad_strata[[i]]) <- names(ssizes) 

gc()

rawbal1_attipdad_strata <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      lapply(
        c(
          "GC-IPW-strata"
        ),
        function(x)
          rawbal1_attipdad_strata[[i]][[j]][[x]]$est   
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
    evidence = "GC-IPW-strata"
  ) 

