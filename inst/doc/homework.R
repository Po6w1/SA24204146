## ----scatter-plot-------------------------------------------------------------
library(ggplot2)
ggplot(mtcars, aes(x=wt, y=mpg)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  labs(title="Scatter Plot of MPG vs Weight", x="Weight (1000 lbs)", y="Miles Per Gallon")

## ----regression-coefficients--------------------------------------------------
model <- lm(mpg ~ wt, data=mtcars)
coef(model)

## ----summary-table------------------------------------------------------------
summary(mtcars$mpg)

## ----boxplot------------------------------------------------------------------
boxplot(mtcars$mpg, main="Boxplot of Miles Per Gallon", ylab="Miles Per Gallon")

## ----correlation-matrix-------------------------------------------------------
cor_matrix <- cor(mtcars)
knitr::kable(cor_matrix, caption="Correlation Matrix of MTCars Dataset")

## ----heatmap------------------------------------------------------------------
heatmap(cor_matrix, symm=TRUE, main="Heatmap of Correlation Matrix")

## -----------------------------------------------------------------------------
# Function to generate random samples from Rayleigh distribution
generate_rayleigh_samples <- function(n, sigma) {
  # Generate uniform random samples
  u <- runif(n)
  # Transform the uniform samples to Rayleigh samples
  rayleigh_samples <- sigma * sqrt(-2 * log(1 - u))
  return(rayleigh_samples)
}

# Parameters
set.seed(123) # For reproducibility
sample_sizes <- c(1000, 5000, 10000) # Different sample sizes
sigmas <- c(1, 2, 3) # Different values for sigma

# Plotting histograms
par(mfrow=c(1, length(sample_sizes)))

for (sigma in sigmas) {
  for (n in sample_sizes) {
    # Generate Rayleigh samples
    samples <- generate_rayleigh_samples(n, sigma)
    
    # Create a histogram
    hist(samples, probability = TRUE, main = paste("Rayleigh(σ =", sigma, ") n =", n),
         xlab = "Value", xlim = c(0, 4 * sigma), breaks = 30)
    
    # Add the theoretical density curve
    x_vals <- seq(0, 4 * sigma, length.out = 100)
    lines(x_vals, (x_vals / sigma^2) * exp(-x_vals^2 / (2 * sigma^2)), col = "blue")
    
    # Calculate and display the mode
    mode_estimate <- max(samples)
    # legend("topright", legend = paste("Mode est.:", round(mode_estimate, 3)), col = "red", lty = 1)
  }
}

# Check theoretical mode value
cat("Theoretical mode for Rayleigh(σ):", sigmas, "\n")


## -----------------------------------------------------------------------------
# Function to generate random samples from a normal location mixture
generate_normal_mixture <- function(n, p1) {
  p2 <- 1 - p1
  mixture_samples <- numeric(n)
  
  for (i in 1:n) {
    if (runif(1) < p1) {
      mixture_samples[i] <- rnorm(1, mean = 0, sd = 1)  # Sample from N(0, 1)
    } else {
      mixture_samples[i] <- rnorm(1, mean = 3, sd = 1)  # Sample from N(3, 1)
    }
  }
  
  return(mixture_samples)
}

# Parameters
set.seed(123) # For reproducibility
sample_size <- 1000
p1_values <- c(0.25, 0.5, 0.75) # Different mixing probabilities

# Set up plotting
# par(mfrow = c(length(p1_values), 1))

for (p1 in p1_values) {
  # Generate samples
  samples <- generate_normal_mixture(sample_size, p1)
  
  # Create histogram
  hist(samples, probability = TRUE, 
       main = paste("Normal Mixture (p1 =", p1, ")"), 
       xlab = "Value", 
       breaks = 30, 
       col = "lightblue", 
       border = "black")
  
  # Overlay density estimation
  dens <- density(samples)
  lines(dens, col = "blue", lwd = 2)
  
  # Add theoretical density components
  x_vals <- seq(-3, 6, length.out = 100)
  theoretical_density <- p1 * dnorm(x_vals, mean = 0, sd = 1) + 
                         (1 - p1) * dnorm(x_vals, mean = 3, sd = 1)
  lines(x_vals, theoretical_density, col = "red", lwd = 2, lty = 2)
  
  legend("topright", legend = c("Kernel Density", "Theoretical Density"), 
         col = c("blue", "red"), lty = 1:2, lwd = 2)
}

# Summary of observations
cat("Observe the plots to analyze bimodality for different values of p1:\n")
cat("Conjecture: The mixture appears bimodal when p1 is around 0.25 to 0.5 and becomes more unimodal as p1 approaches 1.\n")


## -----------------------------------------------------------------------------
# Function to simulate Compound Poisson-Gamma process
simulate_compound_poisson_gamma <- function(lambda, shape, rate, t = 10, n_sim = 10000) {
  # Initialize variables to store results
  X_t_samples <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    # Generate N(t) from Poisson process
    N_t <- rpois(1, lambda * t)
    
    # Generate Y from Gamma distribution for each arrival
    Y_samples <- rgamma(N_t, shape = shape, rate = rate)
    
    # Compute X(t) as the sum of Y_samples
    X_t_samples[sim] <- sum(Y_samples)
  }
  
  # Calculate empirical mean and variance
  empirical_mean <- mean(X_t_samples)
  empirical_variance <- var(X_t_samples)
  
  # Theoretical values
  E_Y1 <- shape / rate                             # Mean of Gamma
  E_Y1_squared <- (shape*(shape + 1)) / (rate^2)  # Variance of Gamma + (E[Y])^2
  
  theoretical_mean <- lambda * t * E_Y1
  theoretical_variance <- lambda * t * E_Y1_squared
  
  # Display results
  cat("Simulated mean X(10):", empirical_mean, "\n")
  cat("Simulated variance X(10):", empirical_variance, "\n")
  cat("Theoretical mean X(10):", theoretical_mean, "\n")
  cat("Theoretical variance X(10):", theoretical_variance, "\n")
}

# Parameters for the simulation
lambda_values <- c(1, 2, 5)          # Different values of lambda
shape_values <- c(2, 3)               # Shape parameter for Gamma
rate_values <- c(1, 1)                # Rate parameter for Gamma

# Simulating for several parameter combinations
for (lambda in lambda_values) {
  for (shape in shape_values) {
    for (rate in rate_values) {
      cat("\nSimulation results for lambda =", lambda, "shape =", shape, "rate =", rate, ":\n")
      simulate_compound_poisson_gamma(lambda, shape, rate)
    }
  }
}


## ----monte-carlo-function-----------------------------------------------------
# Monte Carlo estimate for Beta(3,3) CDF
beta_cdf_mc <- function(x, n_sim = 10000) {
  # Generate random samples from Beta(3,3)
  samples <- rbeta(n_sim, 3, 3)
  # Estimate CDF as proportion of samples <= x
  mean(samples <= x)
}

## ----monte-carlo-estimation---------------------------------------------------
# Values of x
x_values <- seq(0.1, 0.9, 0.1)

# Monte Carlo estimates
mc_estimates <- sapply(x_values, beta_cdf_mc)

# Display the Monte Carlo estimates
mc_estimates

## ----pbeta-comparison---------------------------------------------------------
# Exact values using pbeta
exact_values <- pbeta(x_values, 3, 3)

# Display exact values
exact_values

## ----comparison-table---------------------------------------------------------
# Create a data frame to compare estimates
comparison <- data.frame(
  x = x_values,
  Monte_Carlo_Estimate = mc_estimates,
  Exact_Value = exact_values,
  Difference = abs(mc_estimates - exact_values)
)

# Display the comparison table
comparison

## ----plot-results, fig.width=6, fig.height=4----------------------------------
library(ggplot2)

# Plotting the estimates and exact values
ggplot(comparison, aes(x = x)) +
  geom_point(aes(y = Monte_Carlo_Estimate, color = "Monte Carlo Estimate")) +
  geom_line(aes(y = Monte_Carlo_Estimate, color = "Monte Carlo Estimate")) +
  geom_point(aes(y = Exact_Value, color = "Exact Value")) +
  geom_line(aes(y = Exact_Value, color = "Exact Value")) +
  labs(
    title = "Monte Carlo vs Exact CDF for Beta(3,3)",
    x = "x",
    y = "CDF",
    color = "Method"
  ) +
  theme_minimal()

## ----rayleigh-sample----------------------------------------------------------
# Function to generate Rayleigh samples
rayleigh_sample <- function(sigma, n) {
  U <- runif(n)
  return(sigma * sqrt(-2 * log(U)))
}

## ----antithetic-sample--------------------------------------------------------
# Function to generate Rayleigh samples with antithetic variables
rayleigh_antithetic <- function(sigma, n) {
  U <- runif(n/2)
  X1 <- sigma * sqrt(-2 * log(U))
  X2 <- sigma * sqrt(-2 * log(1 - U))
  return(data.frame(X1 = X1, X2 = X2))
}

## ----variance-estimation------------------------------------------------------
# Set parameters
sigma <- 1  # Rayleigh parameter
n <- 10000  # Number of samples

# Generate independent Rayleigh samples
independent_samples <- data.frame(
  X1 = rayleigh_sample(sigma, n),
  X2 = rayleigh_sample(sigma, n)
)

# Generate Rayleigh samples with antithetic variables
antithetic_samples <- rayleigh_antithetic(sigma, n)

# Calculate means of independent and antithetic samples
mean_independent <- (independent_samples$X1 + independent_samples$X2) / 2
mean_antithetic <- (antithetic_samples$X1 + antithetic_samples$X2) / 2

# Estimate variances
var_independent <- var(mean_independent)
var_antithetic <- var(mean_antithetic)

# Print variances
var_independent
var_antithetic

## ----percent-reduction--------------------------------------------------------
percent_reduction <- (var_independent - var_antithetic) / var_independent * 100
percent_reduction

## ----target-function----------------------------------------------------------
# Target function g(x)
g <- function(x) {
  (x^2 / sqrt(2 * pi)) * exp(-x^2 / 2)
}

## ----importance-function-f1---------------------------------------------------
# Importance function f_1: shifted Gamma distribution
f1 <- function(x, shape = 3, rate = 1) {
  dgamma(x - 1, shape = shape, rate = rate)
}

## ----importance-function-f2---------------------------------------------------
# Importance function f_2: shifted Exponential distribution
f2 <- function(x, rate = 1) {
  dexp(x - 1, rate = rate)
}

## ----sample-f1----------------------------------------------------------------
# Sampling from f_1
set.seed(123)
N <- 10000  # Number of samples
samples_f1 <- rgamma(N, shape = 3, rate = 1) + 1

# Importance sampling estimator using f_1
importance_sampling_f1 <- mean(g(samples_f1) / f1(samples_f1))
importance_sampling_f1

## ----sample-f2----------------------------------------------------------------
# Sampling from f_2
samples_f2 <- rexp(N, rate = 1) + 1

# Importance sampling estimator using f_2
importance_sampling_f2 <- mean(g(samples_f2) / f2(samples_f2))
importance_sampling_f2

## ----variance-comparison------------------------------------------------------
# Variance of the importance sampling estimator for f_1
var_f1 <- var(g(samples_f1) / f1(samples_f1))

# Variance of the importance sampling estimator for f_2
var_f2 <- var(g(samples_f2) / f2(samples_f2))

# Print variances
var_f1
var_f2

## -----------------------------------------------------------------------------
# Quick Sort algorithm implementation
quick_sort <- function(x) {
  if(length(x) <= 1) {
    return(x)
  }
  pivot <- x[1]  # Select the first element as the pivot
  less <- x[x < pivot]  # Elements less than pivot
  equal <- x[x == pivot]  # Elements equal to pivot
  greater <- x[x > pivot]  # Elements greater than pivot
  
  # Recursively sort and combine
  c(quick_sort(less), equal, quick_sort(greater))
}

## ----warning=FALSE, message=FALSE---------------------------------------------
library(microbenchmark)

# Define the sequence of n values
n_values <- c(1e4, 2e4, 4e4, 6e4, 8e4)

# Initialize a list to store the average computation times
avg_times <- numeric(length(n_values))

# Simulate the computation times for each n value
for (i in seq_along(n_values)) {
  n <- n_values[i]
  # Perform 100 simulations and compute average time
  avg_times[i] <- mean(microbenchmark(quick_sort(sample(1:n, n)))$time)
}

# Display the computed average times
avg_times

## -----------------------------------------------------------------------------
# Create a dataframe with n values and corresponding times
data <- data.frame(
  n = n_values,
  t_n = n_values * log(n_values),
  a_n = avg_times
)

# Perform linear regression
model <- lm(a_n ~ t_n, data = data)

# Display the summary of the regression model
summary(model)

## -----------------------------------------------------------------------------
# Create the plot with regression line
ggplot(data, aes(x = t_n, y = a_n)) +
  geom_point(color = "blue", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = "Regression of Computation Time on n log(n)",
    x = expression(t[n] == n * log(n)),
    y = expression(a[n])
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
# Function to calculate skewness
calculate_skewness <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return(m3 / m2^1.5)
}

# Function to generate random normal samples and calculate skewness
generate_skewness_data <- function(m = 10000, n = 50, output_file = "skewness_data.Rds") {
  # Generate skewness statistics
  skstats <- replicate(m, expr = {
    x <- rnorm(n)
    calculate_skewness(x)
  })
  
  # Save the results to file
  saveRDS(skstats, file = output_file)
  
  # Clear memory
  # rm(list = ls())
  
  # Return the file path
  invisible(output_file)
}

# Example of generating and saving skewness data
generate_skewness_data()

## -----------------------------------------------------------------------------
# Function to perform statistical inference on skewness data
perform_inference <- function(input_file = "skewness_data.Rds", output_file = "inference_results.Rds", n = 50, m = 10000) {
  # Load skewness data from file
  skstats <- readRDS(input_file)
  
  # Desired quantiles
  p <- c(0.025, 0.05, 0.95, 0.975)
  
  # Monte Carlo quantiles
  q1 <- quantile(skstats, p)
  
  # Large sample approximation quantiles
  q2 <- qnorm(p, 0, sqrt(6 * (n - 2) / ((n + 1) * (n + 3))))
  q3 <- qnorm(p, 0, sqrt(6 / n))
  
  # Density function for normal approximation
  f <- dnorm(q2, 0, sqrt(6 * (n - 2) / ((n + 1) * (n + 3))))
  
  # Compute standard error of quantile estimates
  v <- p * (1 - p) / (m * f^2)
  se <- sqrt(v)
  
  # Save inference results
  results <- list(
    Percentiles = p,
    MonteCarlo_Quantiles = q1,
    LargeSample_Approx_1 = q2,
    LargeSample_Approx_2 = q3,
    Standard_Error = se
  )
  
  saveRDS(results, file = output_file)
  
  # Clear memory
  # rm(list = ls())
  
  # Return the file path
  invisible(output_file)
}

# Example of performing inference and saving the results
perform_inference()

## -----------------------------------------------------------------------------
# Function to report the results
report_results <- function(input_file = "inference_results.Rds") {
  # Load inference results
  results <- readRDS(input_file)
  
  # Format results into a data frame for reporting
  result_df <- data.frame(
    Percentiles = results$Percentiles,
    MonteCarlo_Quantiles = results$MonteCarlo_Quantiles,
    LargeSample_Approx_1 = results$LargeSample_Approx_1,
    LargeSample_Approx_2 = results$LargeSample_Approx_2,
    Standard_Error = results$Standard_Error
  )
  # Clear memory
  # rm(list = ls())
  
  # Display the results
  knitr::kable(result_df, caption = "Quantile Estimates and Standard Errors")
  
  return(result_df)
}

# Example of reporting the results
report_results()

## -----------------------------------------------------------------------------
# Load necessary libraries
library(MASS)
library(stats) # for cor.test

# Function to generate bivariate normal data
generate_bivariate_normal_data <- function(n = 100, rho = 0.5, m = 1000, output_file = "bivariate_data.Rds") {
  # Generate 'm' bivariate normal samples
  data_list <- replicate(m, {
    mu <- c(0, 0)
    Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    mvrnorm(n, mu, Sigma)
  }, simplify = FALSE)
  
  # Save the generated data
  saveRDS(data_list, file = output_file)
  
  # Clear memory
  # rm(list = ls())
  
  # Return the path to the saved file
  invisible(output_file)
}

# Generate data and save it
generate_bivariate_normal_data()

## -----------------------------------------------------------------------------
# Function to perform correlation tests
perform_correlation_tests <- function(X, Y) {
  pearson_p <- cor.test(X, Y, method = "pearson")$p.value
  spearman_p <- cor.test(X, Y, method = "spearman")$p.value
  kendall_p <- cor.test(X, Y, method = "kendall")$p.value
  
  return(c(pearson_p, spearman_p, kendall_p))
}

# Function to perform inference and compute empirical power
perform_inference <- function(input_file = "bivariate_data.Rds", output_file = "power_results.Rds", alpha = 0.05) {
  # Load the generated data
  data_list <- readRDS(input_file)
  m <- length(data_list)  # Number of simulations
  
  # Perform tests for each dataset
  results <- sapply(data_list, function(data) {
    X <- data[, 1]
    Y <- data[, 2]
    perform_correlation_tests(X, Y)
  })
  
  # Compute empirical power
  pearson_power <- mean(results[1, ] < alpha)
  spearman_power <- mean(results[2, ] < alpha)
  kendall_power <- mean(results[3, ] < alpha)
  
  # Save power results
  power_results_normal <- data.frame(
    Test = c("Pearson", "Spearman", "Kendall"),
    Power = c(pearson_power, spearman_power, kendall_power)
  )
  
  saveRDS(power_results_normal, file = output_file)
  
  # Clear memory
  # rm(list = ls())
  
  # Return the path to the saved file
  invisible(output_file)
}

# Perform inference and save results
perform_inference()

## -----------------------------------------------------------------------------
# Function to report results
report_results <- function(input_file = "power_results.Rds") {
  # Load power results
  power_results_normal <- readRDS(input_file)
  
  # Display the results in a table
  knitr::kable(power_results_normal, caption = "Empirical Power under Bivariate Normal Distribution")
  
  # Clear memory
  # rm(list = ls())
  
  # Return results for further use if needed
  return(power_results_normal)
}

# Report the results
report_results()

## -----------------------------------------------------------------------------
# Function to generate non-normal data
generate_nonnormal_data <- function(n, m, output_file = "data_nonnormal.Rds") {
  results_nonnormal <- replicate(m, {
    X <- runif(n, -1, 1)
    Y <- X^2 + rnorm(n, 0, 0.1)
    cbind(X, Y)
  })
  saveRDS(results_nonnormal, file = output_file)
  # rm(list = ls())  # Clear memory
  invisible(output_file)  # Return the file path
}

# Generate non-normal data and save to file
generate_nonnormal_data(n = 100, m = 1000)

## -----------------------------------------------------------------------------
# Dummy perform_tests function (replace this with your actual test logic)
perform_tests <- function(X, Y) {
  p_pearson <- cor.test(X, Y, method = "pearson")$p.value
  p_spearman <- cor.test(X, Y, method = "spearman")$p.value
  p_kendall <- cor.test(X, Y, method = "kendall")$p.value
  c(p_pearson, p_spearman, p_kendall)
}

# Statistical inference function
perform_inference <- function(input_file = "data_nonnormal.Rds", output_file = "inference_nonnormal.Rds", alpha = 0.05) {
  results_nonnormal <- readRDS(input_file)
  
  m <- dim(results_nonnormal)[3]
  results_test <- replicate(m, {
    data <- results_nonnormal[, , sample(1:m, 1)]  # Select a random sample
    X <- data[, 1]
    Y <- data[, 2]
    perform_tests(X, Y)
  })
  
  # Compute empirical power
  pearson_power_nonnormal <- mean(results_test[1, ] < alpha)
  spearman_power_nonnormal <- mean(results_test[2, ] < alpha)
  kendall_power_nonnormal <- mean(results_test[3, ] < alpha)

  power_results_nonnormal <- data.frame(
    Test = c("Pearson", "Spearman", "Kendall"),
    Power = c(pearson_power_nonnormal, spearman_power_nonnormal, kendall_power_nonnormal)
  )
  
  saveRDS(power_results_nonnormal, file = output_file)
  # rm(list = ls())  # Clear memory
  invisible(output_file)  # Return the file path
}

# Perform inference and save results to file
perform_inference()

## -----------------------------------------------------------------------------
# Result reporting function
report_results <- function(input_file = "inference_nonnormal.Rds") {
  power_results_nonnormal <- readRDS(input_file)
  
  # Display the results in a table
  knitr::kable(power_results_nonnormal, caption = "Empirical Power under Non-normal Alternative")
  # rm(list = ls())  # Clear memory
  return(power_results_nonnormal)  # Return results for further use
}

# Report the results
report_results()

## -----------------------------------------------------------------------------
# Data generation function
generate_data <- function(power1 = 0.651, power2 = 0.676, n_experiments = 10000) {
  output_file <- "data_output.Rds"  # Define the output file inside the function
  data <- list(power1 = power1, power2 = power2, n_experiments = n_experiments)
  saveRDS(data, file = output_file)
  # rm(list = ls())  # Clear memory
  invisible(output_file)  # Return the output file path invisibly
}

# Generate data and save to file
generate_data()


## -----------------------------------------------------------------------------
# Statistical inference function
perform_inference <- function(input_file = "data_output.Rds", output_file = "inference_output.Rds") {
  # Load data
  data <- readRDS(input_file)
  p1 <- data$power1
  p2 <- data$power2
  n <- data$n_experiments
  
  # Compute pooled proportion and standard error
  p_pool <- (p1 * n + p2 * n) / (2 * n)
  se <- sqrt(p_pool * (1 - p_pool) * (1/n + 1/n))
  
  # Z-statistic and p-value
  z_stat <- (p1 - p2) / se
  p_value <- 2 * pnorm(-abs(z_stat))
  
  # Decision
  alpha <- 0.05
  decision <- ifelse(p_value < alpha, "Reject H0: Powers are different", "Fail to Reject H0: No significant difference")
  
  # Save results
  results <- list(Z_statistic = z_stat, P_value = p_value, Decision = decision)
  saveRDS(results, file = output_file)
  # rm(list = ls())  # Clear memory
  invisible(output_file)  # Return the output file path invisibly
}

# Perform inference and save results to file
perform_inference()

## -----------------------------------------------------------------------------
# Result reporting function
report_results <- function(input_file = "inference_output.Rds") {
  # Load results
  results <- readRDS(input_file)
  
  # Report
  knitr::kable(results)
  
  # Clear memory
  # rm(list = ls())
  
  # Optionally, return the results for further use
  return(results)
}

# Report the results
report_results()

## ----simulating-p-values------------------------------------------------------
N <- 1000
n_null <- 950
n_alt <- N - n_null
alpha <- 0.1
m <- 10000 # number of simulation replicates

simulate_pvalues <- function() {
  p_null <- runif(n_null)  # p-values under null hypothesis
  p_alt <- rbeta(n_alt, 0.1, 1)  # p-values under alternative hypothesis
  return(c(p_null, p_alt))
}

## ----applying-corrections-----------------------------------------------------
adjust_pvalues <- function(p_values) {
  bonferroni <- p.adjust(p_values, method = "bonferroni")
  bh <- p.adjust(p_values, method = "BH")
  return(list(bonferroni = bonferroni, bh = bh))
}

## ----calculate-metrics--------------------------------------------------------
calculate_metrics <- function(adjusted_p, true_alt) {
  reject <- adjusted_p < alpha
  true_reject <- true_alt & reject
  false_reject <- !true_alt & reject
  
  fwer <- as.numeric(any(false_reject))
  fdr <- ifelse(sum(reject) > 0, sum(false_reject) / sum(reject), 0)
  tpr <- sum(true_reject) / sum(true_alt)
  
  return(c(fwer, fdr, tpr))
}

## ----run-simulation-----------------------------------------------------------
# Storing results
results_bonf <- matrix(0, nrow = m, ncol = 3)
results_bh <- matrix(0, nrow = m, ncol = 3)

for (i in 1:m) {
  p_values <- simulate_pvalues()
  adjusted_p <- adjust_pvalues(p_values)
  true_alt <- c(rep(FALSE, n_null), rep(TRUE, n_alt))
  results_bonf[i, ] <- calculate_metrics(adjusted_p$bonferroni, true_alt)
  results_bh[i, ] <- calculate_metrics(adjusted_p$bh, true_alt)
}

avg_bonf <- colMeans(results_bonf)
avg_bh <- colMeans(results_bh)

result_table <- matrix(c(avg_bonf, avg_bh), nrow = 3, byrow = FALSE)
colnames(result_table) <- c("Bonferroni Correction", "B-H Correction")
rownames(result_table) <- c("FWER", "FDR", "TPR")
result_table

## ----results-output-----------------------------------------------------------
# Print the table
print(result_table)

## -----------------------------------------------------------------------------
library(boot)
aircond_data <- aircondit$hours

generate_data <- function() {
  return(aircond_data)
}

data_sample <- generate_data()
print(data_sample)  # Display the original data

## -----------------------------------------------------------------------------
mle_hazard_rate <- function(data) {
  lambda_hat <- 1 / mean(data)
  return(lambda_hat)
}

lambda_mle <- mle_hazard_rate(data_sample)
lambda_mle

## -----------------------------------------------------------------------------
bootstrap_stat <- function(data, indices) {
  bootstrap_sample <- data[indices]
  return(1 / mean(bootstrap_sample))
}
R_replicates <- 2000
bootstrap_result <- boot(data_sample, statistic = bootstrap_stat, R = R_replicates)

print(bootstrap_result)

## -----------------------------------------------------------------------------
report_results <- function(bootstrap_obj, original_mle) {
  bias <- mean(bootstrap_obj$t) - original_mle
  se <- sd(bootstrap_obj$t)
  ci <- boot.ci(bootstrap_obj, type = c("norm", "perc", "bca"))
  cat("MLE of Hazard Rate (λ):", original_mle, "\n")
  cat("Bootstrap Bias Estimate:", bias, "\n")
  cat("Bootstrap Standard Error:", se, "\n")
  print(ci)
}

report_results(bootstrap_result, lambda_mle)

## -----------------------------------------------------------------------------
aircond_data <- aircondit$hours

generate_data <- function() {
  return(aircond_data)
}

data_sample <- generate_data()
print(data_sample)  # Display the dataset

## -----------------------------------------------------------------------------
compute_mean <- function(data, indices) {
  sampled_data <- data[indices]  # Resample the data
  return(mean(sampled_data))       # Calculate and return the mean
}

n_replicates <- 2000

bootstrap_result <- boot(
  data = data_sample,               # Input data
  statistic = compute_mean,         # Function for calculating mean
  R = n_replicates                  # Number of bootstrap replicates
)

print(bootstrap_result)

## -----------------------------------------------------------------------------
compute_confidence_intervals <- function(bootstrap_obj) {
  ci_results <- boot.ci(
    boot.out = bootstrap_obj, 
    type = c("norm", "perc", "basic", "bca")  # Different CI methods
  )
  
  cat("Bootstrap Confidence Intervals:\n")
  print(ci_results)
}

compute_confidence_intervals(bootstrap_result)

## -----------------------------------------------------------------------------
report_results <- function(bootstrap_obj) {
  hist(
    bootstrap_obj$t,                # Bootstrap means
    probability = TRUE,              # Density scale
    main = "Distribution of Bootstrap Means",  # Title
    xlab = "Bootstrap Sample Means"  # X-axis label
  )

  points(
    mean(data_sample),              # Original sample mean
    0,                              # Position on the y-axis
    cex = 2,                        # Point size
    pch = 16                        # Point shape (filled circle)
  )
}

report_results(bootstrap_result)

## -----------------------------------------------------------------------------
generate_data_1 <- function() {
  # Load the 'scor' dataset and convert to matrix
  data("scor", package = "bootstrap")
  x <- as.matrix(scor)
  
  # Extract the number of rows (sample size)
  n <- nrow(x)
  
  # Save data to disk for reuse
  save(x, n, file = "data_generation.RData")
  
  cat("Data generation complete. Data saved to 'data_generation.RData'.\n")
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
generate_data_1()   # Generate and save data

## -----------------------------------------------------------------------------
statistical_inference_1 <- function() {
  # Load generated data
  load("data_generation.RData")
  
  # Compute the original statistic (theta_hat)
  lambda <- eigen(cov(x))$values
  theta_hat <- max(lambda / sum(lambda))
  
  # Initialize jackknife estimates
  theta_jack <- numeric(n)
  
  # Perform jackknife resampling
  for (i in 1:n) {
    y <- x[-i, ]  # Leave one observation out
    s <- cov(y)
    lambda <- eigen(s)$values
    theta_jack[i] <- max(lambda / sum(lambda))
  }
  
  # Compute bias and standard error
  bias_jack <- (n - 1) * (mean(theta_jack) - theta_hat)
  se_jack <- sqrt((n - 1) / n * sum((theta_jack - mean(theta_jack))^2))
  
  # Save results to disk
  save(theta_hat, bias_jack, se_jack, file = "statistical_inference.RData")
  
  cat("Statistical inference complete. Results saved to 'statistical_inference.RData'.\n")
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
statistical_inference_1()  # Perform inference and save results

## -----------------------------------------------------------------------------
report_results_1 <- function() {
  # Load inference results from disk
  load("statistical_inference.RData")
  
  # Display the results
  cat("Jackknife Estimation Results:\n")
  cat(sprintf("Theta Hat: %.4f\n", theta_hat))
  cat(sprintf("Jackknife Bias: %.4f\n", bias_jack))
  cat(sprintf("Jackknife SE: %.4f\n", se_jack))
}

## -----------------------------------------------------------------------------
report_results_1()  # Report the saved results

## -----------------------------------------------------------------------------
files_to_delete <- c("data_generation.RData", "statistical_inference.RData")
sapply(files_to_delete, function(file) {
  if (file.exists(file)) {
    file.remove(file)
    cat(file, "has been deleted.\n")
  } else {
    cat(file, "does not exist.\n")
  }
})

## -----------------------------------------------------------------------------
rm(list = ls())  # Clear environment

## -----------------------------------------------------------------------------
generate_data_2 <- function() {
  # Load dataset and convert to vectors
  data("ironslag", package = "DAAG")
  chemical <- ironslag$chemical
  magnetic <- ironslag$magnetic
  
  # Save the data to disk
  save(chemical, magnetic, file = "data_ironslag.RData")
  cat("Data generation complete. Data saved to 'data_ironslag.RData'.\n")
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
generate_data_2()   # Generate and save the data

## -----------------------------------------------------------------------------
statistical_inference_2 <- function() {
  # Load generated data
  load("data_ironslag.RData")
  
  # Create a sequence for predictions
  a <- seq(10, 40, 0.1)

  # Initialize adjusted R-squared values
  Rsq <- numeric(4)

  # Model 1: Linear regression
  L1 <- lm(magnetic ~ chemical)
  Rsq[1] <- summary(L1)$adj.r.squared
  yhat1 <- L1$coef[1] + L1$coef[2] * a

  # Model 2: Quadratic regression
  L2 <- lm(magnetic ~ chemical + I(chemical^2))
  Rsq[2] <- summary(L2)$adj.r.squared
  yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2

  # Model 3: Exponential regression (log-transformed response)
  L3 <- lm(log(magnetic) ~ chemical)
  Rsq[3] <- summary(L3)$adj.r.squared
  logyhat3 <- L3$coef[1] + L3$coef[2] * a
  yhat3 <- exp(logyhat3)

  # Model 4: Cubic regression
  L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))
  Rsq[4] <- summary(L4)$adj.r.squared
  yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 + L4$coef[4] * a^3

  # Perform LOOCV for each model
  n <- length(magnetic)
  e1 <- e2 <- e3 <- e4 <- numeric(n)
  
  for (k in 1:n) {
    # Leave-one-out cross-validation
    y <- magnetic[-k]
    x <- chemical[-k]
    
    # Model 1: Linear
    J1 <- lm(y ~ x)
    e1[k] <- magnetic[k] - (J1$coef[1] + J1$coef[2] * chemical[k])
    
    # Model 2: Quadratic
    J2 <- lm(y ~ x + I(x^2))
    e2[k] <- magnetic[k] - (J2$coef[1] + J2$coef[2] * chemical[k] + 
                            J2$coef[3] * chemical[k]^2)
    
    # Model 3: Exponential
    J3 <- lm(log(y) ~ x)
    e3[k] <- magnetic[k] - exp(J3$coef[1] + J3$coef[2] * chemical[k])
    
    # Model 4: Cubic
    J4 <- lm(y ~ x + I(x^2) + I(x^3))
    e4[k] <- magnetic[k] - (J4$coef[1] + J4$coef[2] * chemical[k] + 
                            J4$coef[3] * chemical[k]^2 + 
                            J4$coef[4] * chemical[k]^3)
  }
  
  # Compute mean squared prediction errors
  mse <- c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

  # Save results to disk
  save(Rsq, mse, file = "inference_results.RData")
  cat("Statistical inference complete. Results saved to 'inference_results.RData'.\n")
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
statistical_inference_2()  # Perform inference and save results

## -----------------------------------------------------------------------------
report_results_2 <- function() {
  # Load inference results from disk
  load("inference_results.RData")
  
  # Display the results
  cat("Adjusted R-squared Values:\n")
  print(Rsq)
  
  cat("\nMean Squared Prediction Errors:\n")
  print(mse)
  
  # Select the best model based on maximum adjusted R^2
  best_model_1 <- which.max(Rsq)
  cat(sprintf("\nThe best model based on maximum adjusted R^2 is Model %d.\n", best_model_1))
  
  # Select the best model based on MSE
  best_model_2 <- which.min(mse)
  cat(sprintf("\nThe model with the lowest prediction error is Model %d.\n", best_model_2))
  cat(sprintf("\nThus the model selected by cross-validation is Model %d.\n", best_model_2))
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
report_results_2()  # Report results

## -----------------------------------------------------------------------------
files_to_delete <- c("data_ironslag.RData", "inference_results.RData")
sapply(files_to_delete, function(file) {
  if (file.exists(file)) {
    file.remove(file)
    cat(file, "has been deleted.\n")
  } else {
    cat(file, "does not exist.\n")
  }
})

## -----------------------------------------------------------------------------
rm(list = ls())  # Clear memory

## -----------------------------------------------------------------------------
# Function: Data Generation
generate_data_3 <- function(feed1, feed2) {
  data(chickwts)     # Load the chickwts dataset
  
  # Extract weights for the two specified feeds
  x <- as.vector(chickwts$weight[chickwts$feed == feed1])
  y <- as.vector(chickwts$weight[chickwts$feed == feed2])
  
  # Save extracted data to disk
  save(x, y, file = "generated_data.RData")
  
  cat(sprintf("Data generated for feeds: %s and %s\n", feed1, feed2))
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
# Function: Statistical Inference (Cramér-von Mises Test)
cvm_permutation_test <- function(R = 500) {
  # Load the generated data
  load("generated_data.RData")
  
  n <- length(x)
  m <- length(y)
  z <- c(x, y)  # Combined data
  N <- n + m    # Total sample size

  # Initialize empirical CDFs
  Fn <- numeric(N)
  Gm <- numeric(N)
  
  # Compute the original test statistic
  for (i in 1:N) {
    Fn[i] <- mean(z[i] <= x)
    Gm[i] <- mean(z[i] <= y)
  }
  cvm0 <- (n * m / N) * sum((Fn - Gm)^2)

  # Perform R permutations
  cvm_permutations <- replicate(R, {
    permuted <- sample(z)  # Random permutation
    X_perm <- permuted[1:n]
    Y_perm <- permuted[(n + 1):N]
    
    # Compute CDFs for permuted data
    for (i in 1:N) {
      Fn[i] <- mean(permuted[i] <= X_perm)
      Gm[i] <- mean(permuted[i] <= Y_perm)
    }
    (n * m / N) * sum((Fn - Gm)^2)
  })

  # Calculate p-value
  p_value <- mean(c(cvm_permutations, cvm0) >= cvm0)

  # Save the inference results
  save(cvm0, p_value, file = "inference_results.RData")
  
  cat("Statistical inference completed.\n")
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
# Function: Result Reporting
report_results_3 <- function() {
  # Load the inference results
  load("inference_results.RData")
  
  # Print the test statistic and p-value
  cat(sprintf("Cramér-von Mises Test Statistic: %.4f\n", cvm0))
  cat(sprintf("Permutation Test P-value: %.4f\n", p_value))
  
  if (p_value < 0.05) {
    cat("The two samples have significantly different distributions (p < 0.05).\n")
  } else {
    cat("No significant difference between the distributions (p >= 0.05).\n")
  }
}


## -----------------------------------------------------------------------------
generate_data_3("soybean", "linseed")
cvm_permutation_test()
report_results_3()

## -----------------------------------------------------------------------------
generate_data_3("soybean", "sunflower")
cvm_permutation_test()
report_results_3()

## -----------------------------------------------------------------------------
generate_data_3("sunflower", "linseed")
cvm_permutation_test()
report_results_3()

## -----------------------------------------------------------------------------
files_to_delete <- c("generated_data.RData", "inference_results.RData")
sapply(files_to_delete, function(file) {
  if (file.exists(file)) {
    file.remove(file)
    cat(file, "has been deleted.\n")
  } else {
    cat(file, "does not exist.\n")
  }
})

## -----------------------------------------------------------------------------
rm(list = ls())  # Clear memory

## -----------------------------------------------------------------------------
# Function: Data Generation
generate_data_4 <- function(n = 30, mu = c(0, 0), rho = 0.5, dist = "normal") {
  # Define the covariance matrix
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  
  # Generate data according to the specified distribution
  if (dist == "normal") {
    data <- mvrnorm(n, mu = mu, Sigma = Sigma)
  } else if (dist == "lognormal") {
    data <- exp(mvrnorm(n, mu = mu, Sigma = Sigma))
  } else {
    stop("Invalid distribution type. Use 'normal' or 'lognormal'.")
  }
  
  # Save the data to disk
  save(data, file = "generated_data.RData")
  cat("Data generated and saved to 'generated_data.RData'.\n")
  
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
generate_data_4(n = 30, rho = 0.5, dist = "normal")

## -----------------------------------------------------------------------------
# Function: Statistical Inference (Spearman Permutation Test)
spear_permutation_test <- function(R = 500) {
  # Load the generated data
  load("generated_data.RData")
  
  x <- data[, 1]
  y <- data[, 2]
  n <- length(x)
  
  # Compute the observed Spearman correlation
  observed_test <- cor.test(x, y, method = "spearman")
  tempp <- observed_test$p.value
  observed_rho <- observed_test$estimate
  
  # Perform R permutations
  permuted_rhos <- replicate(R, {
    y_perm <- sample(y)  # Random permutation of y
    cor.test(x, y_perm, method = "spearman")$estimate
  })
  
  # Combine observed statistic with permutation results
  all_rhos <- c(observed_rho, permuted_rhos)
  
  # Calculate the permutation-based p-value
  perm_p_value <- mean(all_rhos >= observed_rho)
  
  # Save the results
  save(observed_rho, perm_p_value, tempp, file = "inference_results.RData")
  cat("Statistical inference completed and results saved.\n")
  
  rm(list = ls())  # Clear memory
}

## -----------------------------------------------------------------------------
spear_permutation_test()

## -----------------------------------------------------------------------------
# Function: Result Reporting
report_results <- function() {
  # Load the inference results
  load("inference_results.RData")
  
  cat(sprintf("Observed Spearman Correlation: %.4f\n", observed_rho))
  cat(sprintf("Analytical P-value (cor.test): %.4f\n", tempp))
  cat(sprintf("Permutation Test P-value: %.4f\n", perm_p_value))
  
  # Compare the two p-values
  if (abs(perm_p_value - tempp) < 0.05) {
    cat("The p-values from the permutation test and cor.test are similar.\n")
  } else {
    cat("The p-values from the permutation test and cor.test differ significantly.\n")
  }
}

## -----------------------------------------------------------------------------
report_results()

## -----------------------------------------------------------------------------
files_to_delete <- c("generated_data.RData", "inference_results.RData")
sapply(files_to_delete, function(file) {
  if (file.exists(file)) {
    file.remove(file)
    cat(file, "has been deleted.\n")
  } else {
    cat(file, "does not exist.\n")
  }
})

## -----------------------------------------------------------------------------
set.seed(123) # For reproducibility
library(bootstrap)  # For the 'scor' dataset
library(datasets)  # Load base datasets
library(DAAG, warn.conflicts = FALSE)
library(MASS)
library(coda) # For convergence diagnostics (Gelman-Rubin method)

## ----echo=TRUE----------------------------------------------------------------
generate_chains <- function(n_chains = 4, m = 5000, sigma = 3) {
  # Initialize a list to store each chain's samples
  chains <- vector("list", n_chains)

  for (j in 1:n_chains) {
    # Initialize individual chain
    x <- numeric(m)
    x[1] <- rnorm(1, 0, sigma) # Starting value from normal(0, sigma)
    u <- runif(m) # Random uniform numbers for acceptance step

    # Run Metropolis-Hastings algorithm
    for (i in 2:m) {
      xt <- x[i - 1]
      y <- rnorm(1, xt, sigma) # Propose new value

      # Calculate acceptance ratio
      num <- (1 + xt^2) * dnorm(xt, y, sigma)
      den <- (1 + y^2) * dnorm(y, xt, sigma)

      # Accept or reject
      if (u[i] <= num / den) {
        x[i] <- y
      } else {
        x[i] <- xt
      }
    }
    chains[[j]] <- x # Store the chain
  }
  
  return(as.mcmc.list(lapply(chains, function(c) as.mcmc(c))))
  
  rm(list = ls()) # Clear memory
}

## ----echo=TRUE----------------------------------------------------------------
check_convergence <- function(chains) {
  # Apply Gelman-Rubin diagnostic
  r_hat <- gelman.diag(chains)$psrf[1, "Point est."]
  cat("Gelman-Rubin statistic (R-hat):", r_hat, "\n")
  
  return(r_hat < 1.2) # Returns TRUE if converged
  
  rm(list = ls()) # Clear memory
}

## ----echo=TRUE----------------------------------------------------------------
statistical_inference <- function(samples, burn = 1000, 
                                  p = seq(0.1, 0.9, 0.1)) {
  # Extract post-burn-in samples
  xb <- samples[(burn + 1):length(samples)]

  # Compute sample and theoretical quantiles
  Q <- quantile(xb, p)
  theoretical_Q <- qcauchy(p)

  # Display results
  results <- round(rbind(Q, theoretical_Q), 3)
  colnames(results) <- paste0("Decile ", p)
  return(results)
  
  rm(list = ls()) # Clear memory
}

## ----echo=TRUE----------------------------------------------------------------
report_results <- function(samples, burn = 1000) {
  # Extract post-burn-in samples
  xb <- samples[(burn + 1):length(samples)]

  # QQ plot of sample vs theoretical quantiles
  p <- ppoints(100)
  Q <- quantile(xb, p)
  z <- qcauchy(p)

  plot(z, Q, cex = 0.5, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       main = "QQ Plot: Generated Samples vs Theoretical Cauchy")
  abline(0, 1, col = "red")
}

## ----echo=TRUE----------------------------------------------------------------
# Set parameters
n_chains <- 4
m <- 5000
sigma <- 3

# Generate chains
chains <- generate_chains(n_chains, m, sigma)

# Check convergence
converged <- check_convergence(chains)

# If converged, combine chains for further analysis
if (converged) {
  combined_samples <- unlist(as.matrix(chains))
  cat("Chains have converged. Proceeding with statistical inference...\n")
  
  # Perform statistical inference
  decile_comparison <- statistical_inference(combined_samples)
  print("Decile Comparison between Generated Samples and Theoretical Cauchy:")
  print(decile_comparison)

  # Report results with QQ plot
  report_results(combined_samples)
} else {
  cat("Chains have not converged. Please run the chains for more iterations.\n")
}

## -----------------------------------------------------------------------------
generate_rw_laplace <- function(N, x0, sigma) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)  # Uniform samples for acceptance-rejection
  k <- 0  # Counter for rejections
  
  for (i in 2:N) {
    xt <- x[i - 1]
    y <- rnorm(1, mean = xt, sd = sigma)  # Proposal from N(xt, sigma^2)
    
    # Metropolis-Hastings acceptance ratio
    if (u[i] <= exp(abs(xt) - abs(y))) {
      x[i] <- y  # Accept proposal
    } else {
      x[i] <- xt  # Reject proposal
      k <- k + 1  # Increment rejection count
    }
  }
  
  return(list(chain = x, rejections = k))
}

## -----------------------------------------------------------------------------
compute_convergence <- function(chains) {
  mcmc_list <- as.mcmc.list(lapply(chains, mcmc))
  rhat_values <- gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[, 1]
  
  cat("Gelman-Rubin statistics (R-hat):\n", rhat_values, "\n")
  return(rhat_values)
}

## -----------------------------------------------------------------------------
report_results <- function(chains, rejections, N) {
  # Acceptance rates
  acceptance_rates <- 1 - (rejections / N)
  cat("Acceptance Rates:\n", acceptance_rates, "\n")

  # Plot the chains
  # par(mfrow = c(2, 2))
  for (i in seq_along(chains)) {
    plot(chains[[i]], type = "l", main = paste("Chain", i))
  }
  # par(mfrow = c(1, 1))
  
  # Display histograms and compare to Laplace density
  z <- c(-qexp(ppoints(200), 1), qexp(ppoints(200), 1))
  fx <- 0.5 * exp(-abs(z))
  
  # par(mfrow = c(2, 2))
  for (i in seq_along(chains)) {
    hist(chains[[i]], breaks = "Scott", freq = FALSE, ylim = c(0, 0.5), main = paste("Chain", i))
    lines(z, fx, col = "red")
  }
  # par(mfrow = c(1, 1))
}

## -----------------------------------------------------------------------------
set.seed(123)

# Parameters
N <- 5000   # Length of each chain
x0 <- rnorm(1)  # Initial value
sigmas <- c(0.5, 1, 2, 4)  # Proposal standard deviations

# Generate chains for different sigmas
chains <- list()
rejections <- numeric(length(sigmas))

for (i in seq_along(sigmas)) {
  result <- generate_rw_laplace(N, x0, sigmas[i])
  chains[[i]] <- result$chain
  rejections[i] <- result$rejections
}

# Compute Gelman-Rubin statistics
rhat_values <- compute_convergence(chains)

# Run until all chains have R-hat < 1.2
while (any(rhat_values >= 1.2)) {
  cat("Chains not converged. Extending runs...\n")
  
  # Extend the chains
  for (i in seq_along(sigmas)) {
    new_result <- generate_rw_laplace(N, chains[[i]][N], sigmas[i])
    chains[[i]] <- c(chains[[i]], new_result$chain)
    rejections[i] <- rejections[i] + new_result$rejections
  }
  
  N <- N * 2  # Update length
  rhat_values <- compute_convergence(chains)  # Recompute R-hat
}

cat("All chains have converged with R-hat < 1.2.\n")

# Report final results
report_results(chains, rejections, N)

## -----------------------------------------------------------------------------
generate_data <- function(a) {
  # Ensure `a` is a numeric vector
  stopifnot(is.numeric(a))
  
  # Compute the Euclidean norm squared
  norm_sq <- sum(a^2)
  
  # Save norm squared to disk
  saveRDS(norm_sq, file = "norm_squared.rds")
  
  # Return the computed norm squared
  return(norm_sq)
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
compute_sum <- function(K = 60) {
  # Load the norm squared
  norm_sq <- readRDS("norm_squared.rds")
  
  # Initialize required parameters
  k <- 0:K
  d <- length(a)
  
  # Logarithmic calculations for stability with large k and d
  log.ak <- (k + 1) * log(norm_sq)
  log.ck <- lgamma((d + 1) / 2) + lgamma(k + 1.5) - lgamma(k + 1) - 
            k * log(2) - log((2 * k + 1) * (2 * k + 2)) - lgamma(k + d / 2 + 1)
  
  # Compute each term in the sum
  y <- exp(log.ak + log.ck)
  i <- rep(c(1, -1), length = K + 1)
  result <- sqrt(2 / pi) * sum(i * y)
  
  # Save results
  saveRDS(result, file = "sum_result.rds")
  
  # Return the result
  return(result)
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
report_results <- function() {
  # Load the computed result
  sum_result <- readRDS("sum_result.rds")
  
  # Load the norm squared
  norm_sq <- readRDS("norm_squared.rds")
  
  # Report the result
  cat("Computed sum:", sum_result, "\n")
  
  # Approximation for large norm values
  approx_result <- sqrt(norm_sq)
  cat("Approximation (for large norm):", approx_result, "\n")
  
  # Final result (taking the minimum for stability)
  final_result <- min(sum_result, approx_result)
  cat("Final reported result:", final_result)
  
  return(final_result)
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
# Step 1: Data Generation
a <- c(1, 2)
generate_data(a)

# Step 2: Statistical Inference
compute_sum(K = 60)

# Step 3: Result Reporting
report_results()

## -----------------------------------------------------------------------------
generate_data <- function() {
  # Define the range of values for k
  K <- c(4:25, 100, 500, 1000)
  
  # Initial values of a are set to a range for evaluation
  a <- seq(1, sqrt(max(K)) - 0.01, length = 100)
  
  # Save initial data for reuse
  saveRDS(list(K = K, a = a), file = "data_params.rds")
  
  # Return generated data
  return(list(K = K, a = a))
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
solve_for_ck <- function() {
  # Load initial data
  data_params <- readRDS("data_params.rds")
  K <- data_params$K
  a_vals <- data_params$a
  
  # Define helper function to compute CDF difference
  cdf_difference <- function(a, k) {
    c1 <- sqrt(a^2 * (k - 1) / (k - a^2))
    c2 <- sqrt(a^2 * k / (k + 1 - a^2))
    p1 <- pt(c1, df = k - 1, lower.tail = FALSE)
    p2 <- pt(c2, df = k, lower.tail = FALSE)
    p1 - p2
  }
  
  # Initialize storage for results
  a_values <- numeric(length(K))
  pr_values <- numeric(length(K))
  
  # Calculate r values by finding roots where CDF difference is zero
  for (i in seq_along(K)) {
    k <- K[i]
    # Solve for r using uniroot on CDF difference function
    root <- uniroot(cdf_difference, lower = 1, upper = 2, k = k)
    a_values[i] <- root$root
  }
  
  # Calculate c_k based on given formula
  ck <- sqrt(a_values^2 * K / (K + 1 - a_values^2))
  
  # Save c_k results to disk for reporting
  saveRDS(ck, file = "ck_results.rds")
  
  # Save results for reporting
  saveRDS(list(K = K, a = a_values), file = "data_params.rds")
  
  # Return computed values
  return(data.frame(K = K, a = a_values, ck = ck))
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
report_results <- function() {
  # Load computed c_k values
  ck <- readRDS("ck_results.rds")
  
  # Load original data for reference points
  data_params <- readRDS("data_params.rds")
  K <- data_params$K
  a <- data_params$a
  
  # Combine results for display
  results <- cbind(K, a, ck)
  print("Computed results for ck:")
  print(results)
  
  # Return results for further analysis if needed
  return(results)
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
# Step 1: Data Generation
generate_data()

# Step 2: Statistical Inference
solve_for_ck()

# Step 3: Result Reporting
report_results()

## -----------------------------------------------------------------------------
generate_data <- function() {
  # Observed censored data
  Y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
  
  # Censoring threshold
  tau <- 1
  
  # Save initial data to disk
  saveRDS(list(Y = Y, tau = tau), file = "data_params.rds")
  
  # Return data for verification
  return(list(Y = Y, tau = tau))
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
estimate_lambda_em <- function(tol = 1e-6, max_iter = 1000) {
  # Load data parameters
  data_params <- readRDS("data_params.rds")
  Y <- data_params$Y
  tau <- data_params$tau
  
  # Initialize lambda with the mean of the observed data
  lambda <- mean(Y)
  
  # Track convergence
  iter <- 0
  diff <- tol + 1
  
  # EM Algorithm
  while (diff > tol && iter < max_iter) {
    # E-step: Estimate the expected values for censored observations
    uncensored_Y <- Y[Y < tau]
    censored_Y <- Y[Y >= tau]
    n_uncensored <- length(uncensored_Y)
    n_censored <- length(censored_Y)
    
    # Expected total time for censored observations
    E_censored_sum <- n_censored * (tau + 1/lambda)
    
    # M-step: Update lambda
    lambda_new <- (sum(uncensored_Y) + E_censored_sum) / length(Y)
    
    # Check for convergence
    diff <- abs(lambda_new - lambda)
    lambda <- lambda_new
    iter <- iter + 1
  }
  
  # Save final estimate to disk for reporting
  saveRDS(list(lambda = lambda, iter = iter, diff = diff), file = "estimate_results.rds")
  
  # Return the final estimate for verification
  return(list(lambda = lambda, iter = iter, diff = diff))
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
report_results <- function() {
  # Load EM estimate results
  estimate_results <- readRDS("estimate_results.rds")
  lambda_em <- estimate_results$lambda
  iter <- estimate_results$iter
  diff <- estimate_results$diff
  
  # Observed MLE of lambda using uncensored data only
  data_params <- readRDS("data_params.rds")
  Y <- data_params$Y
  tau <- data_params$tau
  uncensored_Y <- Y[Y < tau]
  lambda_mle <- mean(uncensored_Y)
  
  # Display the results
  cat("EM Algorithm Estimation of lambda:\n")
  cat(" - Lambda Estimate (EM):", lambda_em, "\n")
  cat(" - Iterations:", iter, "\n")
  cat(" - Final Convergence Difference:", diff, "\n\n")
  
  cat("Observed Data MLE of lambda (using uncensored values only):", lambda_mle, "\n")
  
  # Plot the estimates for visual comparison
  barplot(c(lambda_em, lambda_mle), beside = TRUE, names.arg = c("EM Estimate", "MLE"),
          col = c("blue", "red"), main = "Comparison of Lambda Estimates",
          ylab = "Lambda Estimate")
  
  # Return results for further use if needed
  return(list(lambda_em = lambda_em, lambda_mle = lambda_mle))
  rm(list = ls()) # Clear memory
}

## -----------------------------------------------------------------------------
# Step 1: Generate data
generate_data()

# Step 2: Estimate lambda using EM algorithm
estimate_lambda_em()

# Step 3: Report and visualize results
report_results()

## -----------------------------------------------------------------------------
library(boot)
library(datasets)
library(bench)
library(fastmatch)

## -----------------------------------------------------------------------------
generate_data <- function() {
  objective_coeffs <- c(4, 2, 9)
  
  constraint_matrix <- rbind(c(2, 1, 1), c(1, -1, 3))
  
  rhs_values <- c(2, 3)
  
  # Return list containing problem data
  list(objective = objective_coeffs, A = constraint_matrix, b = rhs_values)
}

problem_data <- generate_data()
save(problem_data, file = "problem_data.RData")

## -----------------------------------------------------------------------------
apply_simplex <- function() {
  load("problem_data.RData")
  
  solution <- simplex(a = problem_data$objective, A1 = problem_data$A, b1 = problem_data$b, maxi = TRUE)
  
  save(solution, file = "solution.RData")
  solution
}

solution <- apply_simplex()

## -----------------------------------------------------------------------------
report_results <- function() {
  load("solution.RData")
  
  cat("Optimal values of decision variables:\n")
  print(solution$soln) # Decision variables x, y, z
  
  cat("\nMinimum value of the objective function:\n")
  print(solution$value)
}

# Run the result reporting function to display the output
report_results()

## -----------------------------------------------------------------------------
knitr::opts_chunk$set(eval=FALSE, echo = TRUE)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
#  # Version 1: Directly pass the formula list to lapply with lm()
#  models_lapply_v1 <- lapply(formulas, lm, data = mtcars)
#  
#  # Version 2: Use an anonymous function to explicitly pass each formula to lm()
#  models_lapply_v2 <- lapply(formulas, function(formula) lm(formula = formula, data = mtcars))

## -----------------------------------------------------------------------------
#  # Pre-allocate a list to store results from the for loop
#  models_for_loop1 <- vector("list", length(formulas))
#  
#  # Loop through each formula and fit the model
#  for (i in seq_along(formulas)) {
#    models_for_loop1[[i]] <- lm(formulas[[i]], data = mtcars)
#  }

## -----------------------------------------------------------------------------
#  # Display the first model from each approach to compare
#  str(models_lapply_v1[[1]])   # From lapply() version 1
#  str(models_lapply_v2[[1]])   # From lapply() version 2
#  str(models_for_loop1[[1]])    # From for loop

## -----------------------------------------------------------------------------
#  knitr::opts_chunk$set(eval=FALSE, echo = TRUE)
#  bootstraps <- lapply(1:10, function(i) {
#    rows <- sample(1:nrow(mtcars), rep = TRUE)
#    mtcars[rows, ]
#  })

## -----------------------------------------------------------------------------
#  # Generate 10 bootstrap samples of the mtcars dataset
#  set.seed(123) # Set seed for reproducibility
#  bootstraps <- lapply(1:10, function(i) {
#    rows <- sample(1:nrow(mtcars), replace = TRUE) # Sample with replacement
#    mtcars[rows, ] # Return the bootstrap sample
#  })

## -----------------------------------------------------------------------------
#  models_lapply <- lapply(bootstraps, lm, formula = mpg ~ disp)

## -----------------------------------------------------------------------------
#  models_for_loop <- vector("list", length(bootstraps))
#  
#  # Loop through each bootstrap sample and fit the model
#  for (i in seq_along(bootstraps)) {
#    models_for_loop[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
#  }

## -----------------------------------------------------------------------------
#  # Display the structure of the first fitted model from each approach
#  str(models_lapply[[1]])   # Model from lapply() approach
#  str(models_for_loop[[1]]) # Model from for loop approach

## -----------------------------------------------------------------------------
#  knitr::opts_chunk$set(eval=FALSE, echo = TRUE)
#  rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
#  # Define a function to extract R-squared from a model
#  rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
#  # Extract R-squared values from each set of models in Exercise 3
#  r_squared_la1 <- sapply(models_lapply_v1, rsq) # R-squared for models in la1
#  r_squared_la2 <- sapply(models_lapply_v2, rsq) # R-squared for models in la2
#  r_squared_lf1 <- sapply(models_for_loop1, rsq) # R-squared for models in lf1
#  
#  # Display the R-squared values
#  r_squared_la1
#  r_squared_la2
#  r_squared_lf1

## -----------------------------------------------------------------------------
#  # Extract R-squared values from each set of models in Exercise 4
#  r_squared_la <- sapply(models_lapply, rsq) # R-squared for models in la
#  r_squared_lf <- sapply(models_for_loop, rsq) # R-squared for models in lf
#  
#  # Display the R-squared values
#  r_squared_la
#  r_squared_lf

## -----------------------------------------------------------------------------
#  knitr::opts_chunk$set(eval=FALSE, echo = TRUE)
#  trials <- replicate(
#    100,
#    t.test(rpois(10, 10), rpois(7, 10)),
#    simplify = FALSE
#  )

## -----------------------------------------------------------------------------
#  # Set random seed for reproducibility
#  set.seed(123)
#  
#  # Generate 100 trials of t-tests on non-normal data
#  trials <- replicate(
#    100,
#    t.test(rpois(10, 10), rpois(7, 10)),
#    simplify = FALSE
#  )

## -----------------------------------------------------------------------------
#  # Extracting p-values using an anonymous function
#  p_values_anon <- sapply(trials, function(x) x[["p.value"]])
#  
#  # Display the first few p-values
#  head(p_values_anon)

## -----------------------------------------------------------------------------
#  # Extracting p-values without an anonymous function
#  p_values_direct <- sapply(trials, "[[", "p.value")
#  
#  # Display the first few p-values
#  head(p_values_direct)

## -----------------------------------------------------------------------------
#  # Compare the two p-value vectors to confirm they are identical
#  identical(p_values_anon, p_values_direct)

## -----------------------------------------------------------------------------
#  # Sample list of data frames
#  testlist <- list(iris, mtcars, cars)
#  
#  # Inspect the structure of the test list
#  str(testlist)

## -----------------------------------------------------------------------------
#  # Define the custom function using Map() and vapply()
#  lmapply <- function(X, FUN, FUN.VALUE, simplify = FALSE) {
#    # Use Map() to iterate over each element in X, applying vapply() to each
#    out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
#  
#    # If simplify = TRUE, attempt to convert to a matrix or array format
#    if (simplify == TRUE) {
#      return(simplify2array(out))
#    }
#  
#    # Return the result as a list if not simplified
#    out
#  }

## ----warning = FALSE----------------------------------------------------------
#  # Apply the lmapply function to compute column means for each data frame
#  result <- lmapply(testlist, mean, numeric(1))
#  
#  # Display the results
#  print(result)

## ----warning = FALSE----------------------------------------------------------
#  # Apply lmapply with simplify = TRUE
#  result_simplified <- lmapply(testlist, mean, numeric(1), simplify = TRUE)
#  
#  # Display the simplified result
#  print(result_simplified)

## -----------------------------------------------------------------------------
#  chisq.test2 <- function(x, y) {
#    # Combine x and y into a matrix of observed frequencies
#    m <- rbind(x, y)
#  
#    # Calculate row and column sums (margins)
#    margin1 <- rowSums(m)
#    margin2 <- colSums(m)
#    n <- sum(m)
#  
#    # Compute expected frequencies matrix using outer product of margins
#    me <- tcrossprod(margin1, margin2) / n
#  
#    # Calculate chi-squared statistic
#    x_stat <- sum((m - me)^2 / me)
#  
#    # Calculate degrees of freedom
#    df <- (length(margin1) - 1) * (length(margin2) - 1)
#  
#    # Compute p-value
#    p.value <- pchisq(x_stat, df = df, lower.tail = FALSE)
#  
#    # Return results as a list
#    list(x_stat = x_stat, df = df, p.value = p.value)
#  }

## -----------------------------------------------------------------------------
#  # Sample data
#  a <- 21:25
#  b <- seq(21, 29, 2)
#  m <- cbind(a, b)
#  
#  # Results from base chisq.test()
#  base_result <- chisq.test(m)
#  
#  # Results from optimized chisq.test2()
#  opt_result <- chisq.test2(a, b)
#  
#  # Print results for comparison
#  print("Base chisq.test() result:")
#  print(base_result)
#  
#  print("Optimized chisq.test2() result:")
#  print(opt_result)

## -----------------------------------------------------------------------------
#  # Benchmark comparison
#  bench::mark(
#    base = chisq.test(m),
#    optimized = chisq.test2(a, b),
#    check = FALSE
#  )

## -----------------------------------------------------------------------------
#  table2 <- function(a, b) {
#  
#    # Obtain sorted unique values of each vector
#    a_s <- sort(unique(a))
#    b_s <- sort(unique(b))
#  
#    # Calculate dimensions for output table
#    a_l <- length(a_s)
#    b_l <- length(b_s)
#    dims <- c(a_l, b_l)
#    pr <- a_l * b_l
#    dn <- list(a = a_s, b = b_s)
#  
#    # Use fastmatch::fmatch() to quickly map values
#    bin <- fastmatch::fmatch(a, a_s) + a_l * fastmatch::fmatch(b, b_s) - a_l
#  
#    # Use tabulate() to count occurrences
#    y <- tabulate(bin, pr)
#  
#    # Reshape into array format and assign dimension names
#    y <- array(y, dim = dims, dimnames = dn)
#    class(y) <- "table"
#  
#    y
#  }

## -----------------------------------------------------------------------------
#  # Sample input data
#  a <- sample(1:10, 50, replace = TRUE)
#  b <- sample(1:10, 50, replace = TRUE)
#  
#  # Results from standard table()
#  result_table <- table(a, b)
#  
#  # Results from optimized table2()
#  result_table2 <- table2(a, b)
#  
#  # Display comparison
#  identical(result_table, result_table2)

## -----------------------------------------------------------------------------
#  # Larger sample data for benchmarking
#  a <- sample(100, 10000, TRUE)
#  b <- sample(100, 10000, TRUE)
#  
#  # Benchmark
#  bench::mark(
#    base_table = table(a, b),
#    optimized_table2 = table2(a, b),
#    check = FALSE
#  )

## ----cpp_functions------------------------------------------------------------
#  Rcpp::sourceCpp("gibbs_sampler.cpp")

## ----gibbs_example------------------------------------------------------------
#  set.seed(42)
#  # Generate samples using the Gibbs sampler
#  N <- 1000
#  a <- 5
#  b <- 4
#  n <- 20
#  samples_cpp <- gibbs_sampler(N, a, b, n)
#  
#  # Plot the results
#  plot(samples_cpp, main = "Scatter Plot of Gibbs Samples (Rcpp)", xlab = "X", ylab = "Y", col = "blue")

## ----gibbs_r------------------------------------------------------------------
#  gibbs_r <- function(N, a, b, n) {
#    samples <- matrix(0, nrow = N, ncol = 2)
#    x <- 0
#    y <- 0
#  
#    for (i in 1:N) {
#      x <- rbinom(1, n, y)
#      y <- rbeta(1, x + a, n - x + b)
#      samples[i, ] <- c(x, y)
#    }
#    return(samples)
#  }
#  
#  set.seed(42)
#  samples_r <- gibbs_r(N, a, b, n)
#  
#  # Plot the results
#  plot(samples_r, main = "Scatter Plot of Gibbs Samples (Pure R)", xlab = "X", ylab = "Y", col = "red")

## ----qqplot_comparison--------------------------------------------------------
#  qqplot(samples_cpp[, 1], samples_r[, 1], main = "QQ Plot for X (Rcpp vs R)", xlab = "Rcpp", ylab = "R")
#  abline(0, 1, col = "blue")
#  
#  qqplot(samples_cpp[, 2], samples_r[, 2], main = "QQ Plot for Y (Rcpp vs R)", xlab = "Rcpp", ylab = "R")
#  abline(0, 1, col = "blue")

## ----performance_comparison---------------------------------------------------
#  results <- microbenchmark(
#    rcpp = gibbs_sampler(N, a, b, n),
#    pure_r = gibbs_r(N, a, b, n),
#    times = 10
#  )
#  print(results)

## ----uniform_comparison-------------------------------------------------------
#  set.seed(42)
#  cpp_uniform <- cpp_runif(1000, -1, 1)
#  r_uniform <- runif(1000, -1, 1)
#  
#  qqplot(cpp_uniform, r_uniform, main = "QQ Plot for Uniform Samples", xlab = "cpp_runif", ylab = "runif")
#  abline(0, 1, col = "blue")

## ----eval=F-------------------------------------------------------------------
#  d1 <- c(-2.961, 0.478, -0.391, -0.869, -0.460,
#          -0.937, 0.779, -1.409, 0.027, -1.569);
#  d2  <- c(1.608, 1.009,  0.878,  1.600, -0.263,
#           0.680, 2.280,  2.390, 1.793,  8.091, 1.468)

## -----------------------------------------------------------------------------
#  d1 <- c(-2.961, 0.478, -0.391, -0.869, -0.460,
#          -0.937, 0.779, -1.409, 0.027, -1.569)
#  d2 <- c(1.608, 1.009, 0.878, 1.600, -0.263,
#          0.680, 2.280, 2.390, 1.793, 8.091, 1.468)
#  
#  original_mean_diff <- mean(d1) - mean(d2)
#  
#  sample_se <- sqrt(var(d1) / length(d1) + var(d2) / length(d2))
#  
#  R <- 10000
#  bootstrap_diffs <- numeric(R)
#  
#  set.seed(123)
#  for (i in 1:R) {
#    sample_d1 <- sample(d1, length(d1), replace = TRUE)
#    sample_d2 <- sample(d2, length(d2), replace = TRUE)
#    bootstrap_diffs[i] <- mean(sample_d1) - mean(sample_d2)
#  }
#  
#  bootstrap_se <- sd(bootstrap_diffs)
#  
#  cat("Original Mean Difference: ", original_mean_diff, "\n")
#  cat("Sample Standard Error: ", sample_se, "\n")
#  cat("Bootstrap Standard Error: ", bootstrap_se, "\n")

