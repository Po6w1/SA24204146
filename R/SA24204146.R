#' @title Compute and Analyze the Spectrum of a High-Dimensional Random Matrix
#' @name analyze_random_matrix_spectrum
#' @param n Integer. Number of rows in the random matrix.
#' @param p Integer. Number of columns in the random matrix.
#' @param mean Numeric. Mean of the random entries in the matrix. Default is 0.
#' @param sd Numeric. Standard deviation of the random entries. Default is 1.
#' @param bandwidth Numeric. Bandwidth parameter for kernel density estimation. Default is 0.1.
#' @param plot Logical. Whether to plot the results. Default is TRUE.
#' @return A list containing the eigenvalues, theoretical density, and kernel density estimates.
#' @examples
#' library(SA24204146)
#' analyze_random_matrix_spectrum(n = 500, p = 400)
#' @importFrom stats rnorm density
#' @importFrom graphics lines legend
#' @export
analyze_random_matrix_spectrum <- function(n, p, mean = 0, sd = 1, bandwidth = 0.1, plot = TRUE) {
  if (n <= p) stop("Number of rows (n) must be greater than the number of columns (p).")

  # Generate random matrix
  X <- matrix(stats::rnorm(n * p, mean = mean, sd = sd), n, p)

  # Compute the sample covariance matrix
  S <- t(X) %*% X / n

  # Compute eigenvalues of the covariance matrix
  eigenvalues <- eigen(S, only.values = TRUE)$values

  # Theoretical Marcenko-Pastur distribution
  q <- p / n
  lambda_min <- (1 - sqrt(q))^2
  lambda_max <- (1 + sqrt(q))^2
  mp_density <- function(lambda) {
    if (lambda < lambda_min || lambda > lambda_max) return(0)
    sqrt((lambda_max - lambda) * (lambda - lambda_min)) / (2 * pi * q * lambda)
  }

  # Estimate empirical spectral density using kernel density estimation
  density_est <- stats::density(eigenvalues, bw = bandwidth)

  if (plot) {
    # Plot results
    x_seq <- seq(lambda_min, lambda_max, length.out = 1000)
    mp_values <- sapply(x_seq, mp_density)
    plot(density_est, main = "Spectral Density of Random Matrix",
         xlab = "Eigenvalue", ylab = "Density", col = "blue", lwd = 2)
    lines(x_seq, mp_values, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = c("Empirical (KDE)", "Marcenko-Pastur"),
           col = c("blue", "red"), lwd = 2, lty = c(1, 2))
  }

  # Return results
  return(list(
    eigenvalues = eigenvalues,
    mp_density = mp_density,
    kernel_density = density_est
  ))
}

#' @title Monte Carlo Simulation for Diagonal Dominance Proportion
#' @description Simulates the distribution of diagonally dominant row proportions for random matrices.
#' @name monte_carlo_diagonal_dominance
#' @param n Number of rows (and columns) in the random square matrix.
#' @param num_sim Number of Monte Carlo simulations (default is 1000).
#' @param dist Distribution for generating random matrix entries: "normal" (default) or "uniform".
#' @return A ggplot2 histogram and density plot showing the distribution of diagonally dominant row proportions.
#' @examples
#' monte_carlo_diagonal_dominance(n = 5, num_sim = 1000, dist = "normal")
#' monte_carlo_diagonal_dominance(n = 5, num_sim = 1000, dist = "uniform")
#' @import ggplot2
#' @importFrom stats runif
#' @export
utils::globalVariables(c("DominanceRatio", "..density.."))
monte_carlo_diagonal_dominance <- function(n, num_sim = 1000, dist = "normal") {

  dominance_ratios <- numeric(num_sim)  # 存储对角占优行的比例

  for (i in 1:num_sim) {
    # 生成随机矩阵
    if (dist == "normal") {
      mat <- matrix(rnorm(n * n), nrow = n, ncol = n)
    } else if (dist == "uniform") {
      mat <- matrix(runif(n * n, -1, 1), nrow = n, ncol = n)
    } else {
      stop("Unsupported distribution: choose 'normal' or 'uniform'")
    }

    # 使用 Rcpp 函数计算对角占优行的比例
    dominance_ratios[i] <- diagonal_dominance_cpp(mat)
  }

  # 生成对角占优比例分布图
  dominance_df <- data.frame(DominanceRatio = dominance_ratios)

  p <- ggplot(dominance_df, aes(x = DominanceRatio)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(color = "darkblue", size = 1) +
    labs(title = paste("Monte Carlo Distribution of Diagonally Dominant Row Proportions (", dist, ")", sep = ""),
         x = "Proportion of Diagonally Dominant Rows",
         y = "Density") +
    theme_minimal()

  return(p)
}

