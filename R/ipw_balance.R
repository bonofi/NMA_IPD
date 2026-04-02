#' Balance studies using IPW
#' 
#' 


ipw_balance <- function(ipd_network, 
                        model_formula = as.formula(study ~ x + V1 + V2),  # do not balance for X
                        estimand = c("ATE", "ATT"),
                        ref_study = "1",   # for ATT estimation
                        stop_rule = "ks.mean",   # can be a vector
                        n_trees = 3000)
{
  
  estimand <- match.arg(estimand)

  
  # factorize study label
  ipd_network$study <- as.factor(ipd_network$study)

  
  res <- twang::mnps(
    model_formula,
    data = ipd_network |> 
      as.data.frame(),
    estimand = estimand,
    treatATT = ifelse(estimand == "ATE", 
                      NULL, ref_study),
    verbose = FALSE,
    stop.method = stop_rule,
    n.trees = n_trees)
 

  baltable <- twang::bal.table(res)
  
  n <- length(res$psList)
  weights <- lapply(1:n,
                    function(i)
                      data.frame(
                        ps = res$psList[[i]]$ps[res$psList[[i]]$treat == 1, 1],
                        ps_weight = res$psList[[i]]$w[res$psList[[i]]$treat == 1, 1],
                        study = i,
                        usubjid = paste0(i, "-", 1:sum(res$psList[[i]]$treat == 1))
                      )
  )|> 
    dplyr::bind_rows() |> 
    dplyr::mutate(study = as.factor(study))
  
  # print diagnostics
  
  if (!dir.exists("./output/")) 
    dir.create("./output/")
  
  pdf(paste0(
    "./output/inconsistency-", unique(ipd_network$inconsistency),
    "_samplesize-", unique(ipd_network$samplesize), ".pdf"
  ))
  
  plot(res)
  plot(res, plots = 2)
  ggplot(weights,
         aes(study, ps)) + 
    geom_boxplot() + 
    labs(title = "Propensity scores")
  plot(res, plots = 3)
  plot(res, plots = 4)
  plot(res, plots = 5)
  
  grid::grid.newpage()
  # gridExtra::grid.table(
  #   baltable |> 
  #     dplyr::filter(var == "V1")
  # )
  # grid::grid.newpage()
  gridExtra::grid.table(
    baltable |> 
      dplyr::filter(var == "V2")
  )
  dev.off()
  

  # data with weights
  data <- ipd_network |> 
    dplyr::left_join(
      weights,
      by=c("study", "usubjid")
    )

  
  # interpret ATE before-after weighting
  baltable |> 
    filter(var == "V2") |> 
    select(tmt1, mean1, stop.method) |> 
    bind_rows(
      baltable |> 
        filter(var == "V2") |> 
        select(tmt2, mean2, stop.method) |>
        rename(tmt1 = tmt2, mean1 = mean2)
    ) |> distinct() |> 
    rename(
      study = tmt1,
      mean_V2 = mean1
    ) |> 
    mutate(
      mean_V1 = 1 - mean_V2,
      .after = mean_V2
    ) |> 
    arrange(desc(stop.method)) |> 
    left_join(
      data |> 
        distinct(study, trt_name) |> 
        filter(trt_name != "A"),
      by = "study"
    ) |> 
    rowwise() |> 
    mutate(
      weight_eff = case_when(
        trt_name == "B" ~ weighted.mean(
          c(-20,0), c(mean_V2, mean_V1)),
        trt_name == "C" ~ weighted.mean(
          c(-10,0), c(mean_V2, mean_V1)),
        TRUE ~ NA_real_
      )
    ) |> 
    left_join(
      data |> count(study),
      by = "study"
    ) |> 
    group_by(stop.method, trt_name) |> 
    summarise(ATE = weighted.mean(weight_eff, n))
  
  ####  RUN IPWed MODEL ##############
  # weighted model for average effect
  mod <- lm(y~trt_name, data = data, weights = ps_weight)
  
  # (indirect) contrasts
  modcntr <- emmeans::contrast(
    emmeans::emmeans(mod, ~trt_name) ,
    method = "pairwise", 
  )
  modcontr <- modcntr |> 
    as.data.frame() |> 
    dplyr::bind_cols(
      emmeans:::confint.emmGrid(modcntr) |> 
        dplyr::select(dplyr::ends_with("CL"))
    ) |> 
    tibble::add_column(model = NA, .after = Inf) |> 
    dplyr::rename(stat = t.ratio)
  
  
  return(
    list(
      rewres = res,
      weighted_mod = mod,
      ps_weights = data
      est = modcontr |> 
        tibble::add_column(
          evidence = paste0("IPW-", estimand)
        )
    )
  )
  
}




prova <- ipw_balance(
  res1dat |> 
    filter(
      inconsistency == "high",
      samplesize == "small"
    )
)