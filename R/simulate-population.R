
## Function to simulate one population. By construction, the sampling
## is noninformative (inclusion probabilities are independent of
## response)

## One way to introduce dependence is to make the random effects for
## each cluster positively correlated with cluster size (so households
## in larger clusters, which will have lower inclusion probability,
## will also have larger response values)

simulate_survey_population <-
    function(M = 50,          # Number of population clusters
             N_lambda = 100,  # Mean cluster size (Poisson parameter)
             mu = 50,         # Population grand mean
             tau = 5,         # Between-cluster SD
             sigma = 10,      # Within-cluster SD
             rho = 0,         # correlation between random effects and cluster size
             dist = c("normal", "laplace", "lognormal"))
{
    dist <- match.arg(dist)
    
    ## ---------------------------------------------------------
    ## 1. Generate Finite Population
    ## ---------------------------------------------------------
    ## Cluster sizes

    N_i <- rpois(M, lambda = N_lambda)
    ## FIXME: should we ensure N_i >= cluster sample size n?
    ## N_i[N_i < n] <- n + 5

    ## Random effects b_i
    if (rho == 0) {
        b_i <- tau * rnorm(M)
    }
    else {
        ## To make inclusion probabilities correlated with response, a
        ## simple trick is to make the random effects correlated with
        ## cluster size (and hence inclusion probability). But we want
        ## to keep the marginal distributions 

        X <- (N_i - N_lambda) / sqrt(N_lambda) # standardize
        b_i <- tau * (rho * X + sqrt(1 - rho^2) * rnorm(M))
    }


    N_total <- sum(N_i)
    names(N_i) <- as.character(seq_len(M))
    
    ## Generate individual observations
    pop_list <- vector("list", M)
    for (i in seq_len(M)) {
        n_cluster <- N_i[[i]]
        cid <- names(N_i)[[i]] # cluster id
        if (dist == "normal") {
            eps <- rnorm(n_cluster, mean = 0, sd = sigma)
        } else if (dist == "laplace") {
            ## Double exponential: mean = 0, sd = sigma (scale = sigma / sqrt(2))
            u <- runif(n_cluster, -0.5, 0.5)
            eps <- -sign(u) * (sigma / sqrt(2)) * log(1 - 2 * abs(u))
        } else if (dist == "lognormal") {
            ## Skewed/heavy-tailed centered error
            raw <- rlnorm(n_cluster, meanlog = 0, sdlog = 0.5)
            eps <- (raw - exp(0.5^2 / 2)) * (sigma / sqrt(exp(0.5^2) * (exp(0.5^2) - 1)))
        }
        
        pop_list[[i]] <- data.frame(
            cluster_id = cid,
            unit_id = seq_len(n_cluster),
            N_i = n_cluster,
            Y = mu + b_i[i] + eps
        )
    }
    
    do.call(rbind, pop_list)
}


