# -----------------------------
# Set sample size 
# -----------------------------
samples <- c(10, 25, 50, 100, 150, 200)

# -----------------------------
# generate missclas matrix 
# -----------------------------

source("pupulation_model.R")
missclass_matrix_A <-  truncated_misclassification(func_sigma = function(x) 1)
missclass_matrix_B <-  truncated_misclassification(func_sigma = function(x) 4)
missclass_matrix_C <-  truncated_misclassification(func_sigma = function(x) 1, func_mu = function(x) x+4/(x+0.5))
missclass_matrix_D <-  truncated_misclassification(func_sigma = function(x) 4, func_mu = function(x) x+4/(x+0.5))

missclass_all_t <- list(
  A = t(missclass_matrix_A),
  B = t(missclass_matrix_B),
  C = t(missclass_matrix_C),
  D = t(missclass_matrix_D)
)

missclass_all <- list(
  A = missclass_matrix_A,
  B = missclass_matrix_B,
  C = missclass_matrix_C,
  D = missclass_matrix_D
)

# -----------------------------
# simulate samples
# -----------------------------
sim_results_true_dist <- lapply( list(stable_dist), function(dist) {
  lapply(samples, function(n) {
    simulate_counts(dist, size_per_sample = n, n_samples = 100)
  })
})

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

