# this comment

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
                     levels = 1:4,
                     labels = c("A", "B", "C", "D"))
  
  df$sample_size <- factor(df$sample_size,
                           levels = 1:6,
                           labels = c("10", "25", "50", "100", "150", "200"))
  
  df
}
  