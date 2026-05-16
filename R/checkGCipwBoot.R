#' Check GC-IPW bootstrap sample
#' @description
#' check contrasts distribution by scenario, outliers, ....
#' 
#' 

checkGCipwBoot <- function(rawbal1c4)
{
  
  bootsamp <- lapply(
    names(rawbal1c4),
    function(i) lapply(
      names(rawbal1c4[[i]]), 
      function(j)
        rawbal1c4[[i]][[j]][["GC-IPW"]]$rawest |> 
        dplyr::select(contrast, estimate) |> 
        tibble::add_column(
          inconsistency = i,
          samplesize = j)
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
  
  bootsamp |> 
    ggplot2::ggplot(
      aes(x = estimate, fill = samplesize)
    ) + 
    ggplot2::geom_histogram() +
    ggplot2::facet_grid(
      rows = vars(contrast),
      cols = vars(inconsistency)
    )
  
  
}


#checkGCipwBoot(rawbal1c4)