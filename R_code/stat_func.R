source("pupulation_model.R")

# Function to sample from a multinomial distribution
simulate_counts <- function(distr_vector, size_per_sample, n_samples = 1) {
  replicate(
    n_samples,
    as.vector(rmultinom(1, size = size_per_sample, prob = distr_vector))
  )
}
# Total variaton function
tv_distance <- function(p, q) {
  0.5 * sum(abs(p - q))
}

# S simulate convolution while keeping integer counts
simulate_convolved_counts <- function(counts, K) {
  n_bins <- length(counts)
  convolved <- numeric(n_bins)
  
  for (i in seq_len(n_bins)) {
    if (counts[i] > 0) {
      # Sample observed bins for all individuals in true bin i
      convolved <- convolved + rmultinom(1, size = counts[i], prob = K[i, ])
    }
  }
  
  as.vector(convolved)
}

# -----------------------------
# EM algorithm for deconvolution
# -----------------------------

em_deconvolve <- function(counts_convolved, K, n_iter = 100, tol = 1e-4, init = FALSE) {
  
  n_bins <- length(counts_convolved)
  n_total <- sum(counts_convolved)
  
  #  Initialize prior (can also use uniform or your stable_dist)
  if(init){
    p <- stable_dist
  } else {
    p <- rep(1/n_bins, n_bins)
  }
  
  iterations = 1
  for (iter in seq_len(n_iter)) {
    
    # Store old p to check convergence
    p_old <- p
    
    
    # E-step: compute conditoned P(true=i | observed=j)
    posterior <- matrix(0, n_bins, n_bins) # rows=i, cols=j
    for (j in seq_len(n_bins)) {
      denom <- sum(K[, j] * p)    # Total probability of ending in j from any i k is a miscclas-matrix
      if (denom > 0) {
        posterior[, j] <- (K[, j] * p) / denom # P_ij / sum (P_ij)  
      } else {
        posterior[, j] <- rep(1/n_bins, n_bins) # fallback
      }
    }
    
    # M-step: update prior p_i by aggregating contributions over observed bins
    p <- rowSums(sweep(posterior, 2, counts_convolved, `*`)) / n_total
    
    #  Check convergence
    diff <- max(abs(p - p_old))
    rel_diff <- diff / (max(abs(p_old)) + 1e-12)
    # We break if reach convergence or are moving at to slow of a rate
    if (diff < tol || rel_diff < tol) break
    iterations = iterations + 1
  }
  
  # Optional: convert expected counts
  counts_deconvolved <- p * n_total
  
  return(list(p = p, counts = counts_deconvolved, post = posterior, iter = iterations))
}

# Simulate deconvolved counts from observed counts using posterior
simulate_deconvolved_counts <- function(counts_convolved, posterior) {
  
  n_bins <- length(counts_convolved)
  deconvolved <- numeric(n_bins)
  
  for (j in seq_len(n_bins)) {
    if (counts_convolved[j] > 0) {
      # Sample true bins for all individuals in observed bin j
      sampled <- rmultinom(
        n = 1, 
        size = counts_convolved[j], 
        prob = posterior[, j]
      )
      
      deconvolved <- deconvolved + as.vector(sampled)
    }
  }
  
  deconvolved
}

simulate_events <- function(counts, miscclass) {
  true <- rep(seq_along(counts), counts)
  
  obs <- unlist(lapply(seq_along(counts), function(i) {
    sample(seq_along(counts), size = counts[i], replace = TRUE, prob = miscclass[i, ])
  }))
  
  data.frame(true = true, obs = obs, diff = obs - true)
}