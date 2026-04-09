# Function to sample from a discrete distribution
simulate_counts <- function(distr_vector, n) {
  t(rmultinom(n = 1, size = n, prob = distr_vector))
}

# Step 2: simulate convolution while keeping integer counts
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

em_deconvolve <- function(counts_convolved, K, n_iter = 100, tol = 1e-6) {
  
  n_bins <- length(counts_convolved)
  n_total <- sum(counts_convolved)
  
  #  Initialize prior (can also use uniform or your stable_dist)
  p <- rep(1/n_bins, n_bins)
  
  # exp srat
  # bins <- 1:n_bins
  # lambda <- 1 / mean(counts_convolved)
  # p <- dexp(bins, rate = 0.1)
  # p <- p / sum(p)
  
  for (iter in seq_len(n_iter)) {
    
    # Store old p to check convergence
    p_old <- p
    
    
    # E-step: compute posterior P(true=i | observed=j)
    posterior <- matrix(0, n_bins, n_bins) # rows=i, cols=j
    for (j in seq_len(n_bins)) {
      denom <- sum(K[, j] * p)    # Total probability of ending in j from any i 
      if (denom > 0) {
        posterior[, j] <- (K[, j] * p) / denom # P_ij / sum (P_ij)  
      } else {
        posterior[, j] <- rep(1/n_bins, n_bins) # fallback
      }
    }
    
    # M-step: update prior p_i by aggregating contributions over observed bins
    p <- rowSums(sweep(posterior, 2, counts_convolved, `*`)) / n_total
    
    #  Check convergence
    if (max(abs(p - p_old)) < tol) break
  }
  
  # Optional: convert expected counts
  counts_deconvolved <- p * n_total
  
  return(list(p = p, counts = counts_deconvolved, post = posterior))
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

# KL divergence: D_KL(p || q)
kl_divergence <- function(p, q, eps = 1e-12) {
  # Add a small epsilon to avoid log(0)
  p <- p + eps
  q <- q + eps
  
  sum(p * log(p / q))
}