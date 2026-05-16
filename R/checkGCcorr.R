#' check GC cor across all scenarios
#' @description
#' Check correlation diagnostics of GC data using different raw formats ....
#' 


checkGCcorr <- function(rawres, 
                        GConly = FALSE # if GC only input structur is different
){
  
  cordat <- lapply(
    names(rawres),
    function(i) lapply(
      names(rawres[[i]]), 
      function(j)
      {
        if (GConly)
          dat <- lapply(
            names(rawres[[i]][[j]]$raw),
            function(x) 
              rawres[[i]][[j]]$raw[[x]]$is.data.similar$lower.triangular.Rx |> 
              tibble::rownames_to_column("pair") |> 
              dplyr::select(pair, diff) |> 
              tibble::add_column(study = as.integer(x))
          ) 
        else
          dat <- lapply(
            names(rawres[[i]][[j]][["GC-IPW"]]$corr_diagnostics_pseudodat),
            function(x) 
              rawres[[i]][[j]][["GC-IPW"]]$corr_diagnostics_pseudodat[[x]] |> 
              tibble::rownames_to_column("pair") |> 
              dplyr::select(pair, diff) |> 
              tibble::add_column(study = as.integer(x))
          )
        
        dat |>
          dplyr::bind_rows() |> 
          tibble::add_column(
            inconsistency = i,
            samplesize = j)
      }
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
  
  cordat |> 
    ggplot2::ggplot(
      aes(study, diff, colour = samplesize, label = pair)
    ) + 
    ggplot2::geom_point() +
    ggplot2::facet_wrap(vars(inconsistency)) + 
    ggplot2::geom_hline(yintercept = 0) + 
    ggplot2::geom_text(
      data = cordat |> 
        dplyr::filter(abs(diff) >= 0.1),
      hjust = 0, nudge_y = 0.05
    )
   
}





#checkGCcorr(rawGC3)
checkGCcorr(rawbal1c4)