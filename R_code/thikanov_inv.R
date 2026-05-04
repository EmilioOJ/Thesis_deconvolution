# source("simulation_data.R")


# Takes a Matrix of matrixes where col is replika and row is bin and computes TV summary
compute_tv_all <- function(sim_results, stable_dist) {
  
  lapply(seq_along(sim_results), function(i) {
    
    true_dist <- stable_dist
    
    lapply(sim_results[[i]], function(mat) {
      compute_tv_matrix_summary(mat, true_dist) # call compute_tv_matrix_summary
    })
    
  })
}

# Transforms Tv in to summary statistics 
compute_tv_matrix_summary <- function(mat, true_dist) {
  
  tv_vec <- apply(mat, 2, function(counts) {
    est_dist <- counts / sum(counts)
    tv_distance(est_dist, true_dist)
  })
  
  list(
    mean   = round(mean(tv_vec),3),
    median = round(median(tv_vec),3),
    upper_95 = round(quantile(tv_vec, 0.95),3)
  )
}

# creates a df summary
tv_to_df_summary <- function(tv_results) {
  
  out <- list()
  
  idx <- 1
  
  for (i in seq_along(tv_results)) {
    for (j in seq_along(tv_results[[i]])) {
      
      stats <- tv_results[[i]][[j]]
      
      out[[idx]] <- data.frame(
        model = i,
        sample_size = j,
        mean = stats$mean,
        median = stats$median,
        upper_95 = as.numeric(stats$upper_95)
      )
      
      idx <- idx + 1
    }
  }
  df <- do.call(rbind, out)
  
  df$model <- factor(df$model,
                     levels = 1:2,
                     labels = c("σ = 1", "σ = 4"))
  
  df$sample_size <- factor(df$sample_size,
                           levels = 1:5,
                           labels = c("50", "100", "250", "500", "1000"))
  
  df
}



# ------------------------------------------------------------------------------------

# Generate summary for forward simulation


# Outputs a matrix of the projection per simulation n steps in to the future displays it as total population
compute_leslie_projection <- function(mat, leslie_power_list, n_future_proj = 5) {
  
  nr_col <- ncol(mat)
  m <- matrix(NA, nrow = n_future_proj +1, ncol = nr_col)
  
  # iterate over each column
  for (j in seq_len(nr_col) ) {
    
    # First col value 
    obs_pop <- c(sum(mat[,j]))
    
    # where leslie_power_list = list(L,L^2,..,L^n_future_proj)
    for (i in seq_len(n_future_proj) ) {
      obs_pop <- c(obs_pop, sum(leslie_power_list[[i]]%*%mat[,j]))
    }
  m[,j] <- obs_pop
  }
  return(m)
}

# Takes a simulation matrix and sumarises
summary_leslie_proj <- function(mat) {
  
  data.frame(
    mean     = apply(mat, 1, mean),
    median   = apply(mat, 1, median),
    lower_95 = apply(mat, 1, function(x) quantile(x, 0.025)),
    upper_95 = apply(mat, 1, function(x) quantile(x, 0.975))
  )
}


generate_leslie_powers <- function(L, n) {
  
  leslie_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    leslie_list[[i]] <- L %^% i
  }
  
  return(leslie_list)
}

# Compute summary leslie projection nested matrices 

leslie_summary_nested <- function(sim_nested, L, n_future = 5) {
  
  leslie_powers <- generate_leslie_powers(L, n_future)
  
  lapply(sim_nested, function(model_list) {
    
    lapply(model_list, function(mat_list) {
      
      lapply(mat_list, function(mat) {
        
        proj <- compute_leslie_projection(
          mat = mat,
          leslie_power_list = leslie_powers,
          n_future_proj = n_future
        )
        
        summary_leslie_proj(proj)
        
      })
      
    })
    
  })
}

### Median estimator

extract_median_estimator <- function(mat, true_dist) {
  
  tv_vec <- apply(mat, 2, function(counts) {
    est_dist <- counts / sum(counts)
    tv_distance(est_dist, true_dist)
  })
  
  median_index <- which.min(abs(tv_vec - median(tv_vec)))
  
  est_dist <- mat[, median_index]
  est_dist <- est_dist / sum(est_dist)
  
  list(
    estimator = est_dist,
    median_index = median_index,
    tv_value = tv_vec[median_index]
  )
}

get_median_estimator_by_sd_and_n <- function(sd_list, sample_index, true_dist) {
  
  mat <- sd_list[[sample_index]]
  
  extract_median_estimator(mat, true_dist)$estimator
}

get_median_estimator_fast <- function(mat, true_dist) {
  
  # normalize all simulations at once (fast matrix op)
  col_sums <- colSums(mat)
  est_mat <- sweep(mat, 2, col_sums, "/")
  
  # compute TV vector (still loop, but minimal work per step)
  tv_vec <- colSums(abs(est_mat - true_dist)) / 2
  
  median_index <- which.min(abs(tv_vec - median(tv_vec)))
  
  est_mat[, median_index]
}

# Summarise motne carlo per age bin

summarise_mc <- function(mat) {
  data.frame(
    age = seq_len(nrow(mat)),
    median = apply(mat, 1, median),
    lower = apply(mat, 1, quantile, probs = 0.025),
    upper = apply(mat, 1, quantile, probs = 0.975)
  )
}

# Normalizer
summarise_mc <- function(mat, method_name) {
  # normalize each simulation (column)
  mat <- apply(mat, 2, function(x) x / sum(x))
  
  data.frame(
    age = seq_len(nrow(mat)),
    median = apply(mat, 1, median),
    lower = apply(mat, 1, quantile, probs = 0.025),
    upper = apply(mat, 1, quantile, probs = 0.975),
    method = method_name
  )
}