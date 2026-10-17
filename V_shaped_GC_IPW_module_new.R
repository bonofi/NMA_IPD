# run GC method 4 for ATT in level IPD-AGD and AGD
# 
# 

# rawGC3 <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawGC3.rds")
# rawbal1_attipdad <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1_attipdad.rds")
# rawbal1_attad <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1_attad.rds")

######################################## STRATEGY: run GC on optimal settings first, then run IPW on GC data 

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
              method = "11", #  "13",
              SI_k = c(50000, rep(10000, 4)), # only reference study 1 gets the optimal MC setting
              only_SI = TRUE,
              seed = 30697,
              cores = detectCores()-2
            )
            
          }
          
          
        )
    ) 
)


## run IPW IPD-AGD
# must run on subset due to computational resources
set.seed(2607)
index <- sample(1:300, 100)

system.time(
  
  rawbal1_attipdad <- names(rawGC3) |> 
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
            
            list(
              
              "GC-IPW" = gcipdr_ipw_balance(
                ipd_network = res1dat |> 
                  dplyr::filter(inconsistency == i & samplesize == j),
                do_pseudodata = FALSE,
                pseudodata = rawGC3[[i]][[j]]$pseud[index],
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

#####################

allres1b <- res1 |> 
  dplyr::bind_rows(
    res1bbal,
    # res1_attipdad,
    # res1_attad |> 
    #   mutate(
    #     level = "AGD"
    #   ),
    res1_attipdad4,
    res1_attad4
  )



allres1b |> 
  dplyr::filter(evidence2 == "Balanced" & estimand == "ATT") |> 
  dplyr::rowwise() |> 
  dplyr::mutate(Method = paste(evidence, level, sep = "-")) |> 
  dplyr::ungroup() |> 
  ggplot(
    aes(x = samplesize, y = estimate, 
        colour = contrast, group = contrast, shape = Method) 
  ) +
  geom_line(
    data = allres1b |> 
      dplyr::rowwise() |> 
      dplyr::mutate(Method = paste(evidence, level, sep = "-")) |> 
      dplyr::ungroup() |>
      dplyr::filter(Method == "IPW-IPD" & evidence2 == "Balanced" & estimand == "ATT")
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
    title = "ATT estimation in a V-shaped treatment network for increasing inconsistency levels (Nr of studies: 5) after balancing", 
    caption = paste0(
      "Diamonds: IMT for large N"
    )
  ) 


