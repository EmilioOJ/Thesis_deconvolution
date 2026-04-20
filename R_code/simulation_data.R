
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

