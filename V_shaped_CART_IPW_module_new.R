# run CART for ATT in level IPD-AGD and AGD
# 

# rawbal1_attipdad <- readRDS("C:/Users/federico.bonofiglio/Downloads/rawbal1_attipdad.rds")

######################################## STRATEGY: run GC on optimal settings first, then run IPW on GC data 

## run IPW IPD-AGD


rawcart_att <- lapply(
  c("ipd-agd", "agd"),
  \(dl) names(inconsistency) |> 
    purrr::map(
      \(i) names(ssizes) |> 
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
              
              "CART-IPW" = cart_ipw_balance(
                ipd_network = res1dat |> 
                  dplyr::filter(
                    inconsistency == i & samplesize == j),
                modelformula = as.formula(study ~ x + V1 + V2),
                method = "ctree",
                datalevel = dl,
                estimand = "ATT",
                stop_rule = "es.mean",
                cores = detectCores() - 2,
                boot_iter = 10
                
              )
              
            )
          }
          
        )
    ) 
); names(rawcart_att) <- c("ipd-agd", "agd")

for (i in names(rawcart_att)){
  for (j in names(inconsistency)){
    names(rawcart_att[[i]]) <- names(inconsistency)
    names(rawcart_att[[i]][[j]]) <- names(ssizes)
  }
}

gc()


cart_att <- lapply(
  names(rawcart_att),
  \(dl) lapply(
    names(inconsistency),
    function(i) lapply(
      names(ssizes), 
      function(j)
        lapply(
          c(
            "CART-IPW"
          ),
          function(x)
            rawcart_att[[dl]][[i]][[j]][[x]]$est   
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
      evidence = "CART-IPW"
    ) |> 
    bind_rows()
)



################## OLD

system.time(
  
  rawcart_attipdad <- names(inconsistency) |> 
    purrr::map(
      \(i) names(ssizes) |> 
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
              
              "CART-IPW" = cart_ipw_balance(
                ipd_network = res1dat |> 
                  dplyr::filter(inconsistency == i & samplesize == j),
                modelformula = as.formula(study ~ x + V1 + V2),
                method = "ctree",
                estimand = "ATT",
                stop_rule = "es.mean",
                cores = detectCores() - 2,
                boot_iter = 10
                
              )
              
            )
          }
          
        )
    ) 
)

names(rawcart_attipdad) <- names(inconsistency)
for (i in names(inconsistency))
  names(rawcart_attipdad[[i]]) <- names(ssizes) 

saveRDS(rawcart_attipdad, "./results/cartres1.rds")

gc()

cart_attipdad <- lapply(
  names(inconsistency),
  function(i) lapply(
    names(ssizes), 
    function(j)
      lapply(
        c(
          "CART-IPW"
        ),
        function(x)
          rawcart_attipdad[[i]][[j]][[x]]$est   
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
    evidence = "CART-IPW"
  ) 



#####################

allres1b <- res1 |> 
  dplyr::bind_rows(
    res1bbal,
    cart_attipdad  )



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


