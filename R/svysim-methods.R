
## Steps for simulation:

## 1. Generate population (and calculate true mean)
## 2. Generate sample
## 3. Get estimates and conf intervals

## It is not clear whether we should replicate the sampling procedure
## multiple times _within_ a particular generated
## population. Conceivably some methods may perform better for a
## particular sample, but I don't see why knowing that would be
## important. So we will only simulate one sample for every
## randomly generated 'population'.

## This brings up one issue --- the target estimand changes in each
## simulation (as it is the population mean, not the superpopulation
## mean). Thus the sample variance of the estimates have no particular
## meaning as an indicator of performance. We will consider only the
## following:
##
## - bias: average difference (error) between estimand and estimate
## - MSE: average squared error
## - coverage: proportion of cases where estimate +/- 1.96 * SE contains estimand


run_single_simulation <-
    function(M = 50,          # Number of population clusters
             N_lambda = 100,  # Mean cluster size (Poisson parameter)
             mu = 50,         # Population grand mean
             tau = 5,         # Between-cluster SD
             sigma = 10,      # Within-cluster SD
             m = 15,          # Sampled clusters (Stage 1)
             n = 10,          # Sampled units per cluster (Stage 2)
             rho = 0,         # correlation between random effects and cluster size
             dist = c("normal", "laplace", "lognormal"))
{
    dpop <- simulate_survey_population(M = M, N_lambda = N_lambda,
                                       mu = mu, tau = tau,
                                       sigma = sigma,
                                       rho = rho,
                                       dist = dist)
    dsamp <- simulate_survey_sample(dpop, m = m, n = n)
    c(mean(dpop$Y), # true mean
      estimate_survey_mean(dsamp)) # estimates
}


## -----------------------------------------------------------------------------
## Simulation Runner & Performance Metrics (Bias, Variance, MSE)
## -----------------------------------------------------------------------------
svysim <- function(R = 200, ..., cl = NULL)
{
    f <- function(i, ...) run_single_simulation(...) # optionally tryCatch()
    res <-
        if (is.null(cl)) sapply(seq_len(R), f, ..., USE.NAMES = FALSE)
        else parSapply(cl = cl, seq_len(R), f, ..., USE.NAMES = FALSE, chunk.size = 50)
    ## OR
    ## f <- function(i) run_single_simulation(...)
    ## res <- sapply(seq_len(R), f)
    rownames(res) <- c("estimand",
                       "HT", "LME", "PML", "PCL", "Naive",
                       "HT", "LME", "PML", "PCL", "Naive")
    list(estimand  = res[1, ],
         results   = t(res[2:6, ]),
         SE        = t(res[7:11,])) |> structure(class = "svysim")
}


summary.svysim <- function(object, digits = 4, ...)
{
    e          <- with(object, results - estimand)
    ci_ok      <- with(object,
                       results - 1.96 * SE < estimand &
                       results + 1.96 * SE > estimand)
    ci_width   <- with(object, 2 * 1.96 * SE)
    ## Performance evaluation against the true finite population mean Y_bar
    eval_table <-
        data.frame(Bias     = colMeans(e),
                   MSE      = colMeans(e * e),
                   CI_Cov   = colMeans(ci_ok),
                   CI_Width = colMeans(ci_width)) |>
        within(Rel_Eff <- MSE[[1]] / MSE)
    rownames(eval_table) <-
        c("Horvitz-Thompson / Hajek",
          "Unweighted LME (BLUP)",
          "Pseudo-ML (WeMix) BLUP",
          "Pairwise-CL (svylme) BLUP",
          "Naive Weighted Mean")
    round(eval_table, digits)
}

plot.svysim <- function(x, y = NULL, grid = TRUE, ...) {
    d <- cbind(true_mean = x$estimand, as.data.frame(x$results))
    lattice::splom(d, grid = grid, ...)
}




confint.svysim <- function(object,
                           parm = c("bias", "coverage", "width", "mse"),
                           level = 0.95,
                           label = NULL, ...) {
    ## confidence intervals for Bias, Coverage, and CI width, and maybe MSE
    ## do all for now
    e          <- with(object, results - estimand)        [, -5]
    ci_ok      <- with(object,
                       results - 1.96 * SE < estimand &
                       results + 1.96 * SE > estimand)    [, -5]
    ci_width   <- with(object, 2 * 1.96 * SE)             [, -5]
    
    bias <-
        apply(e, 2, t.test, conf.level = level, simplify = FALSE) |>
        sapply("[[", "conf.int") |> t() |>
        as.data.frame() |> setNames(c("bias_lcl", "bias_ucl"))

    cov <-
        colSums(ci_ok) |> sapply(prop.test, n = nrow(ci_ok),
                                 conf.level = level, simplify = FALSE) |>
        sapply("[[", "conf.int") |> t() |>
        as.data.frame() |> setNames(c("cov_lcl", "cov_ucl"))

    ## not really that useful --- closely linked with coverage
    width <-
        apply(ci_width, 2, t.test, conf.level = level, simplify = FALSE) |>
        sapply("[[", "conf.int") |> t() |>
        as.data.frame() |> setNames(c("width_lcl", "width_ucl"))
    
    ## MSE --- what would be a good CI? e^2?

    mse <-
        apply(e^2, 2, t.test, conf.level = level, simplify = FALSE) |>
        sapply("[[", "conf.int") |> t() |>
        as.data.frame() |> setNames(c("MSE_lcl", "MSE_ucl"))
    
    cbind(method = colnames(e),
          label = label %||% "",
          bias, cov, width, mse)

}


