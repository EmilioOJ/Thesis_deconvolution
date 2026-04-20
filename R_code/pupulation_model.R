library(dplyr)
library(tibble)
library(ggplot2)

# computing an actual population model

library(tibble)

L <- matrix(0, 20, 20)
survival <- c(0.5, 0.89, 0.95, 0.99, 0.98, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96)
p_coy <- c(0,0,0,0,0.15,0.35,0.55,0.55,0.55,0.55,0.55,0.7,0.7,0.7,0.7,0.7,0.58,0.58,0.58,0.58)
litter_size <-c(0,0,0,0,2.12,2.25,2.4,2.4,2.4,2.4,2.4,2.6,2.6,2.6,2.6,2.6,2.75,2.75,2.75,2.75) 

F_i <- survival * p_coy * litter_size

length(F_i)

L[1, ] <- F_i

# insert survival in subdiagonal
L[cbind(2:20, 1:19)] <- survival[1:19]  # last survival not used for subdiagonal

# Eigen decomposition
eig <- eigen(L)

# Dominant eigenvalue (largest real part)
lambda_dom <- eig$values[which.max(Re(eig$values))]

# Corresponding right eigenvector
w <- eig$vectors[, which.max(Re(eig$values))]

# Make it real (sometimes tiny imaginary parts appear)
w <- Re(w)

stable_dist <- w / sum(w)

age_labels <- c("0-1", "1-2", "2-3", "3-4", "4-5",
                   "5-6", "6-7", "7-8", "8-9", "9-10",
                   "10-11", "11-12", "12-13", "13-14", "14-15",
                   "15-16", "16-17", "17-18", "18-19", "19+")

#age_labels <- c("pup", as.character(1:19))

# Create dataframe
dist_eigenvector <- data.frame(
  age_labels = age_labels,
  stable_dist = stable_dist
)

# Order factor
dist_eigenvector$age_labels <- factor(
  dist_eigenvector$age_labels,
  levels = age_labels  # <- use the same vector as levels
)

### Ad error

missclass_absorbing <-function(func_sigma= function(x) 1){
  
  error_dist <- numeric(20)
  n_bins <- 20
  missclass_matrix <- matrix(0, nrow = n_bins, ncol = n_bins)
  
  for (i in seq(0.5,19.5,by = 1)){
    
    # heteroskedasticity
    sd <- func_sigma(i)
    # Define intevals
    edges <- c(-Inf, seq(1, 19, by=1), Inf)
    
    # Compute probabilities for each interval
    probs <- pnorm(edges[-1], mean=i, sd=sd) - pnorm(edges[-length(edges)], mean=i, sd=sd)
    
    # Assign to row corresponding to this true age
    row_index <- floor(i + 0.5)  # i = 0.5 → row 1, i = 1.5 → row 2, etc.
    missclass_matrix[row_index, ] <- probs
  }
  return(missclass_matrix)
  
}
## dist with error and truncation

dist_with_error_trunc <-function(X_dist, func_sigma= function(x) 1){
  
  error_dist <- numeric(20)
  
  
  for (i in seq(0.5,19.5,by = 1)){
    
    # heteroskedasticity
    sd <- func_sigma(i)
    # Define intevals
    edges <- c(seq(0, 19, by=1), Inf)
    
    # Compute probabilities for each interval with normalization
    probs <- (pnorm(edges[-1], mean=i, sd=sd) - pnorm(edges[-length(edges)], mean=i, sd=sd))/(1-pnorm(0, mean=i, sd=sd))
    
    # Add contribution of each vector
    error_dist = error_dist + probs*X_dist[i+0.5]
    
  }
  return(error_dist)
  
}
### plot error vs true

plot_overlapping_bars <- function(true_dist, observed_dist,
                                  labels = NULL,
                                  col_true = rgb(1,0,0,0.5),
                                  col_observed = rgb(0,0,1,0.5),
                                  xlab = "Age",
                                  ylab = "Probability",
                                  main = "True vs Observed Distribution") {
  # Number of bars
  n <- length(true_dist)
  
  # Default labels if not provided
  if (is.null(labels)) labels <- 0:(n-1)
  
  # Compute ylim to fit both distributions
  y_max <- max(c(true_dist, observed_dist))
  
  # First barplot: true distribution
  bp <- barplot(height = true_dist,
                names.arg = labels,
                xlab = xlab,
                ylab = ylab,
                main = main,
                col = col_true,
                ylim = c(0, y_max),
                border = NA)
  
  # Second barplot: observed distribution overlaid
  barplot(height = observed_dist,
          col = col_observed,
          border = NA,
          add = TRUE)
  
  # Optional legend
  legend("topright", legend = c("True Distribution", "Observed with Error"),
         fill = c(col_true, col_observed))
  
  invisible(bp)  # return bar midpoints if needed
}


# generate missclasification matrix from a left-truncated N-dist 

truncated_misclassification <- function(func_sigma = function(x) 1, func_mu = function(x) x) {
  
  n_bins <- 20
  missclass_matrix <- matrix(0, nrow = n_bins, ncol = n_bins)
  
  # Define edges of bins: 0–1, 1–2, ..., 19–Inf
  edges <- c(seq(0, n_bins - 1, by = 1), Inf)
  
  for (i in seq(0.5, n_bins - 0.5, by = 1)) {
    
    # heteroskedasticity: standard deviation depends on true age i
    sd <- func_sigma(i)
    # Expected age depend on true age 
    mean_cond_true_age <- func_mu(i)
    
    # Compute raw probabilities for each bin
    probs_raw <- pnorm(edges[-1], mean = mean_cond_true_age, sd = sd) - pnorm(edges[-length(edges)], mean = mean_cond_true_age, sd = sd)
    
    # Normalize to account for truncation at 0 (Y >= 0)
    probs <- probs_raw / (1 - pnorm(0, mean = mean_cond_true_age, sd = sd))
    
    # Assign to row corresponding to this true age
    row_index <- floor(i + 0.5)  # i = 0.5 → row 1, i = 1.5 → row 2, etc.
    missclass_matrix[row_index, ] <- probs
  }
  
  return(missclass_matrix)
}

# Forward modell matrix
forward_model <- function(miss_class, conv_dist) {
  n_bins <- length(conv_dist)
  
  # Column-normalize the misclassification matrix
  f <- miss_class / matrix(colSums(miss_class))
  
  # Redistribute all observed counts
  vec <- f %*% conv_dist
  
  return(vec)
}

# miscclass absolute

missclass_absolute <-function(func_sigma= function(x) 1){
  
  error_dist <- numeric(20)
  n_bins <- 20
  missclass_matrix <- matrix(0, nrow = n_bins, ncol = n_bins)
  
  for (i in seq(0.5,19.5,by = 1)){
    
    # heteroskedasticity
    sd <- func_sigma(i)
    # Define intevals
    edges <- c(-Inf, seq(-18, 19, by=1), Inf)
    
    # Compute probabilities for each interval
    probs_with_neg <- pnorm(edges[-1], mean=i, sd=sd) - pnorm(edges[-length(edges)], mean=i, sd=sd)
    
    # divide in to two vectors
    probs_without_neg <- probs_with_neg[-(1:19)]
    probs_with_neg_only <- probs_with_neg[1:19]
    
    result <- c(probs_without_neg[1], probs_without_neg[2:20] + rev(probs_with_neg_only))
    # Assign to row corresponding to this true age
    row_index <- floor(i + 0.5)  # i = 0.5 → row 1, i = 1.5 → row 2, etc.
    missclass_matrix[row_index, ] <- result
  }
  return(missclass_matrix)
  
}