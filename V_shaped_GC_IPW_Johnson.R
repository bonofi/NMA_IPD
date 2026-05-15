# run GC method 4 for ATT in level IPD-AGD and AGD
# 
# 


system.time(
  
  rawbal1c4 <- split(res1dat, res1dat$inconsistency) |> 
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
                method = "4",
                estimand = "ATT",
                stop_rule = "es.mean",
                boot_iter = 100
              )
              
            )
          }
          
          
        )
    ) 
)

names(rawbal1c4) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1c4[[i]]) <- names(ssizes) 

# rawbal1c4 <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1c4.rds")

res1cbal4 <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      lapply(
        c(
          "GC-IPW"
        ),
        function(x)
          rawbal1c4[[i]][[j]][[x]]$est   
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
    evidence = "GC-IPW4"
  ) 


allres1c <- allres1b |> 
  dplyr::bind_rows(
    res1cbal4
  )


#todo: check if outer level of furrr can speed up module. Check if SI_only is better 

########################################################################
# ---> assuming IPD is NOT available for all studies
# ######################################################################

system.time(
  
  rawbal1d4 <- split(res1dat, res1dat$inconsistency) |> 
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
                method = "4",
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

names(rawbal1d4) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1d4[[i]]) <- names(ssizes) 

# rawbal1d4 <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1d.rds")

res1dbal4 <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      lapply(
        c(
          "GC-IPW"
        ),
        function(x)
          rawbal1d4[[i]][[j]][[x]]$est
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
    evidence = "GC-IPW4"
  )


allres1d <- allres1c4 |> 
  dplyr::bind_rows(
    res1dbal4
  )




######################################## JUST FOCUS ON GC ONCE



system.time(
  
  rawGC4 <- split(res1dat, res1dat$inconsistency) |> 
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
            
              do_gcipdr(
                ipd_network = df2,
                boot_iter = 300,
                method = "4",
                SI_k = 30000,
                only_SI = TRUE,
                seed = 30697,
                cores = detectCores()-2
              )
              
            
          }
          
          
        )
    ) 
)
