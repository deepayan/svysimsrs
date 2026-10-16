
## -----------------------------------------------------------------------------
## Helper Function: royall_blup()
## -----------------------------------------------------------------------------
#' Computes the Royall Finite Population BLUP and its Prediction Variance
#' for a Two-Stage Random Intercept Model:
#'     Y_ij = mu + b_i + eps_ij
#'     b_i ~ iid (0, tau^2),  eps_ij ~ iid (0, sigma^2)
#'
#' @param sample_data A data.frame containing:
#'        - cluster_id: Cluster indicator for sampled units
#'        - Y: Continuous response variable
#' @param N_i_vec A named or ordered numeric vector of all cluster population sizes
#'        in the target universe (length M, names matching cluster_id).
#' @param mu_hat Estimated grand mean / fixed intercept.
#' @param tau2_hat Estimated between-cluster variance component (tau^2).
#' @param sigma2_hat Estimated residual variance component (sigma^2).
#' @param var_mu_hat Optional estimated variance of mu_hat (Var(mu_hat)). 
#'        If NULL, approximated as (tau^2 / m + sigma^2 / (m * n)).
#'
#' @return A list containing:
#'         - est_blup: Point estimate of the finite population mean (Y_bar)
#'         - var_pred: Estimated Mean Squared Error / Prediction Variance of the BLUP
#'         - se_pred:  Standard error of the prediction (sqrt(var_pred))
## -----------------------------------------------------------------------------

royall_blup <-
    function(sample_data,
             N_i_vec,
             mu_hat,
             tau2_hat,
             sigma2_hat,
             var_mu_hat = NULL,
             uhat = NULL)
{
    ## Ensure cluster IDs are treated as character keys
    sample_data$cluster_id <- as.character(sample_data$cluster_id)
    if (is.null(names(N_i_vec)) || !all(sample_data$cluster_id %in% names(N_i_vec))) {
        stop("Names of N_i_vec must match cluster_id")
    }
    all_cluster_ids     <- names(N_i_vec)
    sampled_cluster_ids <- unique(sample_data$cluster_id)
    unsampled_cluster_ids <- setdiff(all_cluster_ids, sampled_cluster_ids)
    
    M <- length(all_cluster_ids)
    m <- length(sampled_cluster_ids)
    N_total <- sum(N_i_vec)
    
    ## Sample sizes and sample sums per sampled cluster
    n_i_vec <- tapply(sample_data$Y, sample_data$cluster_id, length)
    y_bar_i <- tapply(sample_data$Y, sample_data$cluster_id, mean)
    sum_obs <- sum(sample_data$Y)
    
    ## Default variance of fixed effect if not supplied from model output
    if (is.null(var_mu_hat)) {
        mean_n <- mean(n_i_vec)
        var_mu_hat <- (tau2_hat + sigma2_hat / mean_n) / m
    }
    
    ## ---------------------------------------------------------------------------
    ## 1. Point Prediction (BLUP)
    ## ---------------------------------------------------------------------------
    ## Shrinkage factor for each sampled cluster:
    ## gamma_i = tau^2 / (tau^2 + sigma^2 / n_i)
    gamma_vec <- tau2_hat / (tau2_hat + (sigma2_hat / n_i_vec))
    
    ## BLUP for random cluster intercepts: u_hat_i = gamma_i * (y_bar_i - mu_hat)
    u_hat_vec <- gamma_vec * (y_bar_i - mu_hat)
    names(u_hat_vec) <- names(n_i_vec)

    if (!is.null(uhat)) {
        ## FIXME: for lmer, these are different; usually a little, but
        ## sometime a lot --- need to figure out why and which one is correct.
        plot(uhat, u_hat_vec)
        abline(0, 1)
    }
    
    ## Accumulate predictions for unobserved units in SAMPLED clusters:
    ## Unit prediction: Y_hat_ij = mu_hat + u_hat_i
    ## Cluster total unobserved prediction: (N_i - n_i) * (mu_hat + u_hat_i)
    pred_unobs_sampled <- 0
    for (cid in sampled_cluster_ids) {
        n_unobs_i <- N_i_vec[cid] - n_i_vec[cid]
        if (n_unobs_i > 0) {
            pred_unobs_sampled <- pred_unobs_sampled + n_unobs_i * (mu_hat + u_hat_vec[cid])
        }
    }
    
    ## Accumulate predictions for units in UNSAMPLED clusters:
    ## For unsampled clusters, u_hat_k = 0, so Y_hat_kj = mu_hat
    pred_unsampled <- sum(N_i_vec[unsampled_cluster_ids]) * mu_hat
    
    ## Total predicted sum and mean
    total_pred_sum <- sum_obs + pred_unobs_sampled + pred_unsampled
    est_blup <- total_pred_sum / N_total

    ## ---------------------------------------------------------------------------
    ## 2. Prediction Variance / MSE Calculation
    ## ---------------------------------------------------------------------------
    ## The finite population prediction error is:
    ##     e_pred = Y_bar_hat - Y_bar = (1 / N_total) * (T_hat - T)
    ##
    ## Since observed units cancel out exactly (sum_{i in s} Y_i - sum_{i in s} Y_i = 0),
    ## T_hat - T decomposes into unobserved units across two independent sets:
    ##
    ##   T_hat - T = sum_{i in s_cl} (N_i - n_i) * [ (mu_hat + u_hat_i) - (mu + b_i) - e_bar_unobs,i ]
    ##             + sum_{k in uns_cl} N_k * [ mu_hat - (mu + b_k) - e_bar_k ]
    ##
    ## Under Royall (1970) and BHF (1988), the variance components of Var(T_hat - T) are:
    ##
    ## (A) Unsampled Clusters (k in uns_cl):
    ##     - True cluster mean has variance: Var(b_k) = tau^2
    ##     - Unobserved individuals have residual variance: N_k * sigma^2
    ##     - Contribution = sum_{k} [ N_k^2 * tau^2 + N_k * sigma^2 ]
    ##
    ## (B) Sampled Clusters (i in s_cl):
    ##     - Uncertainty in predicting b_i via BLUP u_hat_i:
    ##       Var(u_hat_i - b_i) = tau^2 * (1 - gamma_i) = tau^2 * sigma^2 / (n_i * tau^2 + sigma^2)
    ##     - Unobserved individuals have residual variance: (N_i - n_i) * sigma^2
    ##     - Contribution = sum_{i} [ (N_i - n_i)^2 * tau^2 * (1 - gamma_i) + (N_i - n_i) * sigma^2 ]
    ##
    ## (C) Parameter Estimation Uncertainty for mu_hat:
    ##     - Total unobserved population size: N_unobs_total = sum_{k uns} N_k + sum_{i s} (N_i - n_i) * (1 - gamma_i)
    ##     - Contribution = [ N_unobs_effective ]^2 * Var(mu_hat)
    ## ---------------------------------------------------------------------------
    
    ## Component A: Unsampled clusters
    N_uns <- N_i_vec[unsampled_cluster_ids]
    var_comp_unsampled <- sum(N_uns^2 * tau2_hat + N_uns * sigma2_hat)
    
    ## Component B: Sampled clusters
    var_comp_sampled <- 0
    N_unobs_eff_sampled <- 0
    for (cid in sampled_cluster_ids) {
        n_unobs_i <- N_i_vec[cid] - n_i_vec[cid]
        if (n_unobs_i > 0) {
            g_i <- gamma_vec[cid]
            ## Variance of prediction error of random intercept (tau^2 * (1 - gamma_i))
            var_cluster_effect <- tau2_hat * (1 - g_i)
            var_comp_sampled <- var_comp_sampled + (n_unobs_i^2 * var_cluster_effect + n_unobs_i * sigma2_hat)
            N_unobs_eff_sampled <- N_unobs_eff_sampled + n_unobs_i * (1 - g_i)
        }
    }
    
    ## Component C: Uncertainty due to estimating mu
    N_unobs_effective <- sum(N_uns) + N_unobs_eff_sampled
    var_comp_mu <- (N_unobs_effective^2) * var_mu_hat
    
    ## Total prediction variance of finite population total T
    var_total_pred <- var_comp_unsampled + var_comp_sampled + var_comp_mu
    
    ## Scaled to finite population mean: Var(Y_bar_hat - Y_bar) = Var(T_hat - T) / N_total^2
    var_ybar_blup <- var_total_pred / (N_total^2)
    se_ybar_blup  <- sqrt(var_ybar_blup)
    
    return(c(
        est = as.numeric(est_blup),
        SE  = as.numeric(se_ybar_blup)
    ))
}


## Separate functions for different methods for easy
## plug-and-play. Each function takes argument 'dsample' as returned
## by simulate_survey_sample() and returns estimated mean and SE.


emean_HT <- function(dsample, design = NULL)
{
    ## N_i <- attr(dsample, "N_i") # if needed for fpc
    des <-
        design %||% svydesign(id = ~ cluster_id + unit_id,
                              ## weights = ~ weight, # equivalent for svymean()
                              weights = ~ w1 + w2, # needed for PCL
                              ## fpc = ~rep(M, nrow(dsample)) + N_i,
                              data = dsample)
    fm <- svymean(~ Y, design = des)
    c(est = coef(fm)[[1]],
      SE = SE(fm)[[1]])
}


emean_LME <- function(dsample)
{
    ## FIXME: Needed for Royall BLUP, but can we compute aggregate from weights?
    N_i <- attr(dsample, "N_i")
    fm <- suppressMessages(lmer(Y ~ 1 + (1 | cluster_id), data = dsample, REML = TRUE))
    mu_hat <- fixef(fm)[[1]]
    tau2_hat <- VarCorr(fm)[[1]][[1]] # variance
    sigma_hat <- sigma(fm)       # std. dev
    ## Royall finite population predictor:
    royall_blup(sample_data = dsample,
                N_i_vec = N_i,
                mu_hat = mu_hat,
                tau2_hat = tau2_hat,
                sigma2_hat = sigma_hat^2, # square std. dev
                var_mu_hat = NULL)
}

emean_WLME <- function(dsample)
{
    N_i <- attr(dsample, "N_i")
    dsample$weight[] <- dsample$weight * (nrow(dsample) / sum(dsample$weight))
    fm <- suppressMessages(lmer(Y ~ 1 + (1 | cluster_id), data = dsample,
                                weights = weight, REML = TRUE))
    mu_hat <- fixef(fm)[[1]]
    tau2_hat <- VarCorr(fm)[[1]][[1]] # variance
    sigma_hat <- sigma(fm)       # std. dev
    ## Royall finite population predictor:
    royall_blup(sample_data = dsample,
                N_i_vec = N_i,
                mu_hat = mu_hat,
                tau2_hat = tau2_hat,
                sigma2_hat = sigma_hat^2, # square std. dev
                var_mu_hat = NULL)
}

emean_PML <- function(dsample, ...)
{
    N_i <- attr(dsample, "N_i")
    fm <- mix(Y ~ 1 + (1 | cluster_id),
              data = dsample,
              weights = c("w2", "w1"), # c(level-1, level-2) weights
              ...,
              verbose = FALSE)
    mu_hat <- coef(fm)[[1]]
    ## BLUP random effects under weighted model
    sigma2_hat <- fm$vars[[2]] # ["Residual"]
    tau2_hat <- fm$vars[[1]] # ["cluster_id.(Intercept)"]
    ## Royall finite population predictor:
    royall_blup(sample_data = dsample,
                N_i_vec = N_i,
                mu_hat = mu_hat,
                tau2_hat = tau2_hat,
                sigma2_hat = sigma2_hat,
                var_mu_hat = NULL) 
}


emean_PCL <- function(dsample, design = NULL)
{
    N_i <- attr(dsample, "N_i")
    des <-
        design %||% svydesign(id = ~ cluster_id + unit_id,
                              ## weights = ~ weight, # equivalent for svymean()
                              weights = ~ w1 + w2, # needed for PCL
                              ## fpc = ~rep(M, nrow(dsample)) + N_i,
                              data = dsample)
    fm <- svy2lme(Y ~ 1 + (1 | cluster_id), design = des)
    mu_hat <- coef(fm)[[1]]
    sigma2 <- fm$s2[[1]]
    tau2 <- sigma2 * diag(fm$L)
    ## Royall finite population predictor:
    royall_blup(sample_data = dsample,
                N_i_vec = N_i,
                mu_hat = mu_hat,
                tau2_hat = tau2,
                sigma2_hat = sigma2,
                var_mu_hat = NULL)
}

emean_Naive <- function(dsample, design = NULL)
{
    fm <- lm(Y ~ 1, dsample, weights = weight)
    c(est = coef(fm)[[1]],
      SE = SE(fm)[[1]])
}


estimate_survey_mean <- function(dsample)
{
    ## svydesign --- needed for HT and PCL
    ## N_i <- attr(dsample, "N_i") # if needed for fpc
    des <- svydesign(
        id = ~ cluster_id + unit_id,
        ## weights = ~ weight, # equivalent for svymean()
        weights = ~ w1 + w2, # needed for PCL
        ## fpc = ~rep(M, nrow(dsample)) + N_i,
        data = dsample)

    ## ---------------------------------------------------------
    ## Estimators
    ## 1: Horvitz-Thompson / Hajek Estimator
    ## 2: Unweighted Mixed Model (BHF / Royall BLUP)
    ## 3: Weighted Mixed Model / Pseudo-ML (BLUP)
    ## 4: Pairwise Composite Likelihood
    ## 5: Naive (sample mean)
    ## ---------------------------------------------------------

    rbind(HT = emean_HT(dsample, des),
          LME = emean_LME(dsample),
          ## PML = emean_WLME(dsample), # fast
          ## OR
          PML = emean_PML(dsample), # results differ from WLME but slow
          PCL = emean_PCL(dsample, des),
          Naive = emean_Naive(dsample))

}

