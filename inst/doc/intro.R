## -----------------------------------------------------------------------------
library(SA24204146)

# Analyze the spectrum of a random matrix with 500 rows and 400 columns
result <- analyze_random_matrix_spectrum(n = 500, p = 400, bandwidth = 0.1)

## -----------------------------------------------------------------------------
library(SA24204146)

# Example 1: Normal distribution
set.seed(123)
p1 <- monte_carlo_diagonal_dominance(n = 5, num_sim = 1000, dist = "normal")
print(p1)

# Example 2: Uniform distribution
set.seed(123)
p2 <- monte_carlo_diagonal_dominance(n = 5, num_sim = 1000, dist = "uniform")
print(p2)

## -----------------------------------------------------------------------------
# Define the range of the theoretical Marcenko-Pastur distribution
q <- 400 / 500  # p / n
lambda_min <- (1 - sqrt(q))^2
lambda_max <- (1 + sqrt(q))^2

# Compute the Marcenko-Pastur density over the theoretical range
x_seq <- seq(lambda_min, lambda_max, length.out = 1000)
mp_values <- sapply(x_seq, result$mp_density)

# Plot the results
plot(result$kernel_density, main = "Spectral Density of Random Matrix",
     xlab = "Eigenvalue", ylab = "Density", col = "blue", lwd = 2)
lines(x_seq, mp_values, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Empirical (KDE)", "Marcenko-Pastur"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))



## -----------------------------------------------------------------------------
library(SA24204146)

# Generate and visualize the distribution
set.seed(123)
p <- monte_carlo_diagonal_dominance(n = 5, num_sim = 1000, dist = "normal")
print(p)

