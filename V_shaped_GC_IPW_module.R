# run GC method 4 for ATT in level IPD-AGD and AGD
# 
# 

# rawGC3 <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawGC3.rds")
# rawGC4 <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawGC4.rds")
# rawbal1c4 <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1c4.rds")
# rawbal1d4 <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1d4.rds")
# rawbal1c <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1c.rds")
# rawbal1d <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1d.rds")
# rawbal1_attipdad <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1_attipdad.rds")
# rawbal1_attad <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1_attad.rds")
# rawbal1_attipdad4 <- readRDS("C:/Users/federico.bonofiglio.CYTEL/Downloads/rawbal1_attipdad4.rds")
# rawbal1_attad4 <- readRDS("C:/Users/federico.bonofiglio.CYTEL/Downloads/rawbal1_attad4.rds")
# 

######################################## STRATEGY: run GC on optimal settings first, then run IPW on GC data 

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
            
            gc()
            
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


## 


system.time(
  
  rawGC3 <- split(res1dat, res1dat$inconsistency) |> 
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
            
            do_gcipdr(
              ipd_network = df2,
              boot_iter = 300,
              method = "3",
              SI_k = 30000,
              only_SI = TRUE,
              seed = 30697,
              cores = detectCores()-2
            )
            
          }
          
          
        )
    ) 
)


#### RUN GC-IPW IPD-AD (reference study nr 1 is in IPD format, all other are AD)

rawGC <- rawGC3
# for IPD-AD application, substitute GC study 1 with "original" reference study 1
for (i in names(rawGC))
  for (j in names(rawGC[[i]]))
    for (b in 1:length(rawGC[[i]][[j]]$pseud))
      rawGC[[i]][[j]]$pseud[[b]][which(
        rawGC[[i]][[j]]$pseud[[b]]$study == "1"), ] <- res1dat |> 
  dplyr::filter(study == "1" & inconsistency == i & samplesize == j) |> 
  dplyr::select(
    any_of(colnames(rawGC[[i]][[j]]$pseud[[b]]))
  )


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


#### RUN GC-IPW AD (reference study is nr 1) all studies have AD format


system.time(
  
  rawbal1_attad <- names(rawGC3) |> 
    purrr::map(
      \(i) names(rawGC3[[i]]) |> 
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
                ipd_network = rawGC3[[i]][[j]]$pseud[sample(1:300, 100)],
                do_pseudodata = FALSE,
                modelformula = as.formula(study ~ x + V1 + V2),
                estimand = "ATT",
                datalevel = "agd",
                stop_rule = "es.mean",
                cores = detectCores() - 2
                
              )
              
            )
          }
          
        )
    ) 
)

names(rawbal1_attad) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1_attad[[i]]) <- names(ssizes) 

gc()


res1_attad <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      lapply(
        c(
          "GC-IPW"
        ),
        function(x)
          rawbal1_attad[[i]][[j]][[x]]$est   
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


checkGCipwBoot(rawbal1_attad)

########################## IPW on GC4 !!!!!! ############################
########################## 

#### RUN GC-IPW IPD-AD (reference study nr 1 is in IPD format, all other are AD)

rawGC <- rawGC4
# for IPD-AD application, substitute GC study 1 with "original" reference study 1
for (i in names(rawGC))
  for (j in names(rawGC[[i]]))
    for (b in 1:length(rawGC[[i]][[j]]$pseud))
      rawGC[[i]][[j]]$pseud[[b]] <- res1dat |> 
  dplyr::filter(study == "1" & inconsistency == i & samplesize == j) |> 
  dplyr::select(
    any_of(colnames(rawGC[[i]][[j]]$pseud[[b]]))
  ) |> 
  bind_rows(
    rawGC[[i]][[j]]$pseud[[b]][which(
      rawGC[[i]][[j]]$pseud[[b]]$study != "1"), ]
  )


## run IPW

system.time(
  
  rawbal1_attipdad4 <- names(rawGC) |> 
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

names(rawbal1_attipdad4) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1_attipdad4[[i]]) <- names(ssizes) 

gc()

res1_attipdad4 <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      lapply(
        c(
          "GC-IPW"
        ),
        function(x)
          rawbal1_attipdad4[[i]][[j]][[x]]$est   
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
    evidence = "GC-IPW-4"
  ) 


#### RUN GC-IPW AD (reference study is nr 1) all studies have AD format


system.time(
  
  rawbal1_attad4 <- names(rawGC4) |> 
    purrr::map(
      \(i) names(rawGC4[[i]]) |> 
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
                ipd_network = rawGC4[[i]][[j]]$pseud[sample(1:300, 100)],
                do_pseudodata = FALSE,
                modelformula = as.formula(study ~ x + V1 + V2),
                estimand = "ATT",
                datalevel = "agd",
                stop_rule = "es.mean",
                cores = detectCores() - 2
                
              )
              
            )
          }
          
        )
    ) 
)

names(rawbal1_attad4) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawbal1_attad4[[i]]) <- names(ssizes) 

gc()


res1_attad4 <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      lapply(
        c(
          "GC-IPW"
        ),
        function(x)
          rawbal1_attad4[[i]][[j]][[x]]$est   
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
    evidence = "GC-IPW-4"
  ) 


checkGCipwBoot(rawbal1_attad)

#####################

allres1b <- res1 |> 
  dplyr::bind_rows(
    res1bbal,
    res1_attipdad,
    res1_attad |> 
      mutate(
        level = "AGD"
      ),
    res1_attipdad4,
    res1_attad4
  )


