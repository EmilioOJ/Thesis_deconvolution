# -----------------------------
# Set sample size 
# -----------------------------
samples <- c(50, 100, 250, 500, 1000)

# -----------------------------
# generate missclas matrix 
# -----------------------------

source("pupulation_model.R") # uses misclass models
missclass_matrix_A <-  truncated_misclassification(func_sigma = function(x) 1)
missclass_matrix_B <-  truncated_misclassification(func_sigma = function(x) 4)

missclass_all_t <- list(
  sd1 = t(missclass_matrix_A),
  sd4 = t(missclass_matrix_B)
)

missclass_all <- list(
  sd1 = missclass_matrix_A,
  sd4 = missclass_matrix_B
)

# -----------------------------
# simulate samples
# -----------------------------
sim_results_true_dist <- lapply( list(stable_dist), function(dist) {
  lapply(samples, function(n) {
    simulate_counts(dist, size_per_sample = n, n_samples = 1000)
  })
})
# Normalized estimate 


normalized_truedist_sim <- lapply(sim_results_true_dist, function(model_list) {
  
  lapply(model_list, function(mat) {
    
    stopifnot(is.matrix(mat))
    
    if (anyNA(mat)) stop("NA detected in matrix")
    
    col_sums <- colSums(mat)
    
    if (any(col_sums <= 0)) {
      stop("Invalid column sums (<=0) detected")
    }
    
    sweep(mat, 2, col_sums, "/")
  })
  
})

saveRDS(normalized_truedist_sim, "../data/simulations/normalized_truedist_sim.rds")


# -----------------------------
# Add error to sample 
# -----------------------------

convolved_sim <- lapply(missclass_all, function(K) {
  
  lapply(sim_results_true_dist, function(sim_list) {
    
    lapply(sim_list, function(mat) {
      
      apply(mat, 2, simulate_convolved_counts, K = K)
      
    })
    
  })
  
})

# normalized convolved sim 

normalized_convolved_sim <- lapply(convolved_sim, function(model_list) {
  
  lapply(model_list, function(sample_list) {
    
    lapply(sample_list, function(mat) {
      
      stopifnot(is.matrix(mat))
      
      col_sums <- colSums(mat)
      
      # avoid division by zero
      col_sums[col_sums == 0] <- 1
      
      sweep(mat, 2, col_sums, `/`)
    })
    
  })
  
})

saveRDS(normalized_convolved_sim, "../data/simulations/normalized_convolved_sim.rds")

# -----------------------------
# Deconvolution with ML
# -----------------------------

deconvolved_ml_p_samples_list <- Map(function(model_list, K) {
  
  lapply(model_list, function(scenario_list) {
    
    lapply(scenario_list, function(mat) {
      
      stopifnot(is.matrix(mat))  # IMPORTANT DEBUG SAFETY
      
      apply(mat, 2, function(col) {
        em_deconvolve(col, K = K, init=TRUE, n_iter =1000)
      })
      
    })
    
  })
  
}, convolved_sim, missclass_all)

deconvolved_ml_p_samples_list_false_init <- Map(function(model_list, K) {
  
  lapply(model_list, function(scenario_list) {
    
    lapply(scenario_list, function(mat) {
      
      stopifnot(is.matrix(mat))  # IMPORTANT DEBUG SAFETY
      
      apply(mat, 2, function(col) {
        em_deconvolve(col, K = K, init=FALSE, n_iter =1000)
      })
      
    })
    
  })
  
}, convolved_sim, missclass_all)

# Result
iter_only_u <- lapply(deconvolved_ml_p_samples_list_false_init, function(sd) {
  lapply(sd, function(scenario) {
    lapply(scenario, function(samples) {
      lapply(samples, `[[`, "iter")
    })
  })
})

deconvolved_ml_p_samples_u <- lapply(deconvolved_ml_p_samples_list_false_init, function(sd) {
  lapply(sd, function(scenario) {
    lapply(scenario, function(samples) {
      
      do.call(cbind, lapply(samples, `[[`, "p"))
      
    })
  })
})

saveRDS(iter_only_u, "../data/simulations/iter_no_stop.rds")
saveRDS(deconvolved_ml_p_samples_u, "../data/simulations/pi_samples_ml_u_no_stop.rds")

# -----------------------------
# Deconvolution with MoM
# -----------------------------


lambdas <- list(
  sd1 = 1e-1,
  sd4 = 1e-2
)

deconvolved_mom_p_samples_i <- Map(
  
  function(model_list, K, name) {
    
    lambda <- lambdas[[name]]
    
    lapply(model_list, function(scenario_list) {
      
      lapply(scenario_list, function(mat) {
        
        apply(mat, 2, function(col) {
          
          # normalize counts to proportions
          if (sum(col) > 0) {
            col <- col / sum(col)
          }
          
          # Step 1: regularized inversion (Tikhonov)
          p_hat <- solve(
            t(K) %*% K + lambda * diag(ncol(K)),
            t(K) %*% col
          )
          
          # Step 2: remove negatives
          p_hat <- pmax(p_hat, 0)
          
          # Step 3: normalize
          s <- sum(p_hat)
          if (s > 0) {
            p_hat <- p_hat / s
          }
          
          as.vector(p_hat)
        })
        
      })
      
    })
    
  },
  
  convolved_sim,
  missclass_all,
  names(missclass_all)
)

saveRDS(deconvolved_mom_p_samples_i, "../data/simulations/deconvolved_mom_p_samples_i.rds")