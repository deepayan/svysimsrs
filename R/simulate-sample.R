
## Function to run one simulation iteration

simulate_survey_sample <-
    function(dpopulation,     # data frame giving population
             m = 15,          # Sampled clusters (Stage 1)
             n = 10           # Sampled units per cluster (Stage 2)
             )
{
    ## ---------------------------------------------------------
    ## 1. Generate Finite Population
    ## ---------------------------------------------------------
    ## Extract necessary population characteristics: number of
    ## clusters M and cluster sizes N_i
    dclusters <- unique(dpopulation[c("cluster_id", "N_i")])
    M <- nrow(dclusters)
    N_i <- dclusters$N_i
    names(N_i) <- dclusters$cluster_id
    
    ## ---------------------------------------------------------
    ## 2. Two-Stage Sampling Scheme
    ## ---------------------------------------------------------
    ## Stage 1: Sample m clusters out of M using SRSWOR
    sampled_clusters <- sample(names(N_i), size = m, replace = FALSE)
    pi_1 <- m / M
    
    ## Stage 2: Sample n units from each selected cluster using SRSWOR
    sample_list <- vector("list", m)
    for (idx in seq_along(sampled_clusters)) {
        cid <- sampled_clusters[idx]
        cluster_sub <- dpopulation[dpopulation$cluster_id == cid, ]
        sampled_rows <- sample(seq_len(nrow(cluster_sub)), size = n, replace = FALSE)
        sample_list[[idx]] <- cluster_sub[sampled_rows, ]
    }
    sample_data <- do.call(rbind, sample_list)
    
    ## Survey Weights
    ## Stage 1 weight: w_1i = 1 / pi_1 = M / m
    ## Stage 2 weight: w_2ij = 1 / pi_2|i = N_i / n
    ## Overall inclusion probability: pi_ij = (m / M) * (n / N_i)
    sample_data$w1 <- M / m
    sample_data$w2 <- sample_data$N_i / n
    sample_data$weight <- sample_data$w1 * sample_data$w2
    sample_data$cluster_id <- as.factor(sample_data$cluster_id)
    sample_data |> structure(N_i = N_i,
                             Y_bar_true = mean(dpopulation$Y))
}
