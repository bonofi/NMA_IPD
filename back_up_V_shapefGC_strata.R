
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
              
              # out$level_1$pseud[[b]] |> 
              # dplyr::select(
              #   !dplyr::starts_with("V")
              # ) |> 
              # tibble::add_column(V1 = 1, V2 = 0) |>  # replace with bang-bang and diagonal 0-1 V matrix ?
              # dplyr::bind_rows(
              #   out$level_2$pseud[[b]] |> 
              #     dplyr::select(
              #       !dplyr::starts_with("V")
              #     ) |> 
              #     tibble::add_column(V1 = 0, V2 = 1)
              # ) |> 
              # dplyr::arrange(study)
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
