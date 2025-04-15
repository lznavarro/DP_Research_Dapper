library(MASS)
library(stats)
library(pracma)


rHadamard <- function(d, t=3) {
  
  H <- hadamard(d) / sqrt(d)
  
  random_sign_flip <- function(H) {
    flip <- sample(c(-1, 1), d, replace = TRUE)
    D <- diag(flip)
    return(H %*% D)
  }
  
  HD <- random_sign_flip(H)
  for (i in 2:t) {
    HD <- random_sign_flip(H) %*% HD
  }
  #HD <- t(HD) #suspicious
  return(HD)
}

#' Compute Log Density for Noisy Quantile Binary Search (QBS)
dQBS_sx <- function(svect, counts, m, u, p, T = 10, l) {
  left <- l
  right <- u
  lden_sum <- 0
  
  for (i in 1:T) {
    mid <- (left + right) / 2
    count <- counts[i]
    
    # Add log-density contribution for the current iteration
    lden_sum <- lden_sum + dnorm(svect[i], mean = count, sd = sqrt(T / (2 * p)), log = TRUE)
    
    # Update the search interval based on the noisy count
    if (svect[i] <= m) {
      left <- mid
    } else {
      right <- mid
    }
  }
  
  return(lden_sum)
}

#' Compute Log Density for Clipped Mean with Summary Statistics
d_clipped_mean_sx <- function(noisy_mean, svect, C, sx, n, d, u, p, l = 0, threshold = NULL) {
  p1 <- p * 0.25
  p2 <- p * 0.75
  
  if (!is.null(threshold)) {
    lden_sum <- 0
  } else {
    # Perform QBS to compute log density for counts
    T <- 10
    m <- as.integer(n - 2 * sqrt(d / (2 * p2)) - sqrt(T / (2 * p1)))
    lden_sum <- dQBS_sx(svect, sx$counts_clipping, m, u, p1, T = T, l = 0)
  }
  
  # Compute the log density of the noisy mean
  noise_std <- 2 * C / sqrt(2 * p2) / n
  dnorm_LD_sum <- sum(dnorm(noisy_mean, mean = sx$clipped_mean, sd = noise_std, log = TRUE))
  
  # Add the log density of the noisy mean to the total
  lden_sum <- lden_sum + dnorm_LD_sum
  
  return(lden_sum)
}

#' This function applies a randomized Hadamard transform, clips data points exceeding
#' the threshold, and computes the mean of the clipped data.
S_clipped_mean <- function(x, sdp) {
  
  # Apply Hadamard transformation to the input data
  x_hat <- t(sdp$HD %*% t(x))
  
  # Normalize the transformed data by subtracting the offsets
  Cs_vector <- as.numeric(unlist(sdp$Cs))
  x_shifted <- x_hat - matrix(Cs_vector, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  
  y_list <- list()
  
  for (i in 1:nrow(x_shifted)) {
    xi_norm <- sqrt(sum(x_shifted[i, ]^2))  # Compute the norm of the current data point
    
    # Clip the data point if it exceeds the threshold
    if (xi_norm <= sdp$C_clipped) {
      yi <- x_shifted[i, ]
    } else {
      yi <- (sdp$C_clipped / xi_norm) * x_shifted[i, ]
    }
    
    y_list[[i]] <- yi  
  }
  
  # Compute the mean of the clipped data
  S <- Reduce(`+`, y_list)
  S <- S / sdp$n
  
  return(S)
}


#' Perform Quantile Binary Search (QBS) for sufficient statistics
#' This function uses a binary search approach, guided by noisy observations (`svects`),
#' to iteratively refine the interval and compute counts.
s_QBS <- function(x, svects, u, m, l) {
  T <- length(svects)  # Number of binary search iterations
  
  left <- l
  right <- u
  counts <- rep(0, T) 
  
  for (j in 1:T) {
    mid <- (left + right) / 2  # Compute midpoint of the interval
    count <- length(x[x <= mid])  # Count elements less than or equal to `mid`
    counts[j] <- count  # Store the count for this iteration
    
    # Adjust the search interval based on noisy count
    if (svects[j] <= m) {
      left <- mid
    } else {
      right <- mid
    }
  }
  
  return(counts)
}

#' This function transforms the input data using a randomized Hadamard matrix, 
#' computes clipped statistics, and estimates sufficient statistics for use in privacy-preserving computations.
s_instance_optimal <- function(x, sdp) {
  
  HD <- sdp$HD
  Cs <- sdp$Cs
  svects <- sdp$svects
  u <- sdp$u
  n <- sdp$n
  p <- sdp$p
  T <- sdp$T
  svect_clipped <- sdp$svect_clipped
  counts_list <- list()
  d <- sdp$d

  p <- p * 0.75
  p1 <- p * 0.25
  p2 <- p * 0.75
  m <- as.integer(n - 2 * sqrt(d / (2 * p2)) - sqrt(T / (2 * p1)))
  
  # Transform data using the randomized Hadamard matrix
  x_hat <- t(HD %*% t(x))
  
  # Compute counts for each dimension using Quantile Binary Search
  for (j in 1:d) {
    x_dim <- x_hat[, j, drop = FALSE]
    counts_medians <- s_QBS(x_dim, svects[[j]], u, n / 2, l = -u)
    counts_list[[j]] <- counts_medians
  }
  
  # Shift the data based on precomputed clipping thresholds
  if (nrow(x) == 1) {
    x_shifted <- x_hat - matrix(as.numeric(unlist(Cs)), nrow = 1, ncol = d, byrow = TRUE)
  } else {
    x_shifted <- x_hat - matrix(as.numeric(unlist(Cs)), nrow = nrow(x_hat), ncol = d, byrow = TRUE)
  }
  
  # Compute norms of shifted data
  x_norm <- apply(x_shifted, 1, function(row) sqrt(sum(row^2)))
  
  # Compute counts after applying clipping
  counts_clipping <- s_QBS(x_norm, svect_clipped, u, m, l = 0)
  
  # Compute the clipped mean
  S_clipped <- S_clipped_mean(x, sdp)
  
  s_instance_result <- list(
    counts_list = counts_list,
    counts_clipping = counts_clipping,  
    clipped_mean = S_clipped            
  )
  
  return(s_instance_result)
}

rQBS <- function(x, m, u, p, T=10, l) {
  left <- l
  right <- u
  svect <- c()
  #print(paste("rQBS M:", m))
  counts_list <- numeric(T)
  for (i in 1:T) {
    mid <- (left + right) / 2
    count <- length(x[x <= mid])
    counts_list[i] <- count
    count_noisy <- count + rnorm(1, 0, sqrt(T / (2 * p)))
    svect <- c(svect, count_noisy)
    if (count_noisy <= m) {
      left <- mid
    } else {
      right <- mid
    }
    #print(paste("rQBS Mid:", mid, "Left:", left, "Right:", right))
  }
  #print(svect)
  #print(paste("Counts from rQBS:", paste(counts_list, collapse = ", ")))
  return(list(svect = svect, final_value = (left + right) / 2))
}

r_clipped_mean <- function(x, n, d, u, p, l=0, threshold=NULL) {
  p1 <- p * 0.25
  p2 <- p * 0.75
  x_norm <- apply(x, 1, function(row) sqrt(sum(row^2)))
  if (!is.null(threshold)) {
    C <- threshold
    svect <- NULL
  } else {
    T <- 10
    m <- as.integer(n - 2 * sqrt(d / (2 * p2)) - sqrt(T / (2 * p1)))
    #print(paste("R n:", n, "d:", d, "p1:", p1, "p2:", p2))
    QBS_out <- rQBS(x_norm, m, u, p1, T=T, l=0)
    C <- QBS_out[[2]]
    svect <- QBS_out[[1]]
  }
  #cat("C (clipping threshold):", C, "\n")
  #cat("Proportion of x_norm greater than C:", sum(x_norm > C) / length(x_norm), "\n")
  x_clipped <- list()
  for (i in 1:length(x_norm)) {
    xi_norm <- x_norm[i]
    scale <- min(C / xi_norm, 1.0)
    x_clipped[[i]] <- scale * x[i, ]
  }
  mean <- colMeans(do.call(rbind, x_clipped))
  #cat('The mean is:', mean, '\n')
  
  noise_std <- 2 * C / sqrt(2 * p2) / n
  noisy_mean <- mean + rnorm(d, 0, noise_std)
  
  #noisy_mean <- noisy_mean + mean(x)
  
  return(list(noisy_mean, svect, C))
}

rInstanceOptimal <- function(x, n, d, u, p, T=10, l=1) {
  p1 <- p * 0.25 / d
  p2 <- p * 0.75
  
  HD <- rHadamard(d, t=3)
  HD_inv <- ginv(HD)
  #HD_inv <- t(HD)
  #print(dim(HD))
  #print(dim(x))
  
  x_hat <- t(HD %*% t(x))
  
  svects <- list()
  Cs <- list()
  noisy_mean <- list()
  #cat('Colmeans of x_hat', colMeans(x_hat), '\n')
  
  for (i in 1:d) {
    x_dim <- x_hat[, i, drop = FALSE]
    QBS_out <- rQBS(x_dim, n/2, u, p1, T=T, l=-u)
    svects[[i]] <- QBS_out$svect
    Cs[[i]] <- QBS_out$final_value
  }
  
  #cat('Cs: from rIO', as.numeric(Cs), '\n')
  
  x_shifted <- x_hat - matrix(as.numeric(unlist(Cs)), nrow = n, ncol = d, byrow = TRUE)
  
  clipped_mean_out <- r_clipped_mean(x_shifted, n, d, u, p2, l=0, threshold=NULL)
  
  Cs_matrix <- matrix(as.numeric(unlist(Cs)), nrow = d, ncol = 1) 
  
  clipped_mean_vec <- matrix(clipped_mean_out[[1]], nrow = d, ncol = 1) 
  
  noisy_mean <- as.vector(HD_inv %*% (clipped_mean_vec + Cs_matrix))
  
  return(list(
    noisy_mean = noisy_mean,
    svect_clipped = clipped_mean_out[[2]],
    C_clipped = clipped_mean_out[[3]],
    svects = svects,
    Cs = Cs,
    x = x,
    HD = HD,
    n = n,
    d = d,
    u = u,
    p = p,
    T = T,
    l = l
  ))
}

#' Compute sufficient statistics for a given observation
st_f <- function(xi, sdp, i) {
  
  # Ensure xi is a matrix; convert if necessary
  if (is.null(nrow(xi))) {
    x <- matrix(xi, nrow = 1, ncol = length(xi))
  } else {
    x <- xi
  }
  
  # Compute sufficient statistics using the s_instance_optimal function
  summary_result <- s_instance_optimal(x, sdp)
  
  return(summary_result)
}


#' Calculate the Log-Density of the Privacy-Preserving Mechanism
priv_f <- function(sdp, sx) {
  noisy_mean <- randomVars$noisy_mean
  svect_clipped <- randomVars$svect_clipped
  C_clipped <- randomVars$C_clipped
  svects <- randomVars$svects
  Cs <- randomVars$Cs
  n <- randomVars$n
  HD <- randomVars$HD
  d <- randomVars$d
  u <- randomVars$u
  p <- randomVars$p
  T <- randomVars$T
  l <- randomVars$l
  
  counts_list <- sx$counts_list
  counts_clipping <- sx$counts_clipping
  counts_mean <- sx$counts_mean
  
  # Compute privacy budget allocation
  p1 <- p * 0.25 / d
  p2 <- p * 0.75
  m <- as.integer(n - 2 * sqrt(d / (2 * p2)) - sqrt(T / (2 * p1)))
  
  total <- 0
  
  # Compute the log-density for each dimension using QBS
  for (i in 1:d) {
    svect <- svects[[i]]  
    total <- total + dQBS_sx(svect, counts_list[[i]], m, u, p1, T=T, l=-u)
  }
  
  # Reconstruct the original noisy mean from the Hadamard transformation
  noisy_mean_original <- (HD %*% noisy_mean) - as.numeric(Cs)
  
  # Add log-density from the clipped mean
  CM_lden <- d_clipped_mean_sx(noisy_mean_original, svect_clipped, C_clipped, sx, n,
                               1, u, p2, l=l, threshold=NULL)
  total <- total + CM_lden
  
  return(total)
}

#' Merge Two Summary Objects with Specified Operation
merge_sx <- function(sx1, sx2, add = TRUE) {

  op <- if (add) `+` else `-`

  result_counts_list <- list()
  result_counts_clipping <- numeric(length(sx1$counts_clipping))
  result_clipped_mean <- numeric(length(sx1$clipped_mean))
  
  # Merge counts list element-wise
  for (i in seq_along(sx1$counts_list)) {
    result_counts_list[[i]] <- op(sx1$counts_list[[i]], sx2$counts_list[[i]])
  }
  
  # Apply the operation to clipped counts and clipped means
  result_counts_clipping <- op(sx1$counts_clipping, sx2$counts_clipping)
  result_clipped_mean <- op(sx1$clipped_mean, sx2$clipped_mean)
  
  return(list(
    counts_list = result_counts_list,
    counts_clipping = result_counts_clipping,
    clipped_mean = result_clipped_mean
  ))
}

n <- 200
d <- 4
u <- 20 # keep fixed
p <- 1
mu <- 4
sa <- 2
x <- matrix(rnorm(n * d, mean = mu, sd = sa), ncol = d)

randomVars <- rInstanceOptimal(x, n, d, u, p)

sdp <- list(
  HD = randomVars$HD,
  Cs = randomVars$Cs,
  C_clipped = randomVars$C_clipped,
  svect_clipped = randomVars$svect_clipped,
  noisy_mean = randomVars$noisy_mean,
  n = n,
  svects = randomVars$svects,
  u = randomVars$u,
  p = randomVars$p,
  T = randomVars$T,
  x = randomVars$x,
  d = randomVars$d
)