#TPC_code
# TPC_Code
Thermal Performance Curve, fitting models
################################Example_with_rTPC and nls.multstart
# === STEP 1: Install and load packages ===
# install.packages(c("rTPC", "nls.multstart", "terra", "ggplot2", "tidyverse"))
library(rTPC)
library(nls.multstart)
library(terra)
library(ggplot2)
library(tidyverse)

# === STEP 2: Create example TPC data ===
# This simulates experimental data measuring performance across temperatures
# Replace this with your actual data!

set.seed(123)

# Example: Mosquito development rate data
tpc_data <- data.frame(
 temp = c(10, 12, 15, 18, 20, 22, 25, 27, 30, 32, 35, 37, 40),
 performance = c(0.05, 0.12, 0.25, 0.45, 0.65, 0.85, 0.95, 0.92, 0.75, 0.55, 0.30, 0.10, 0.01)
)

# Add some measurement error (realistic data has noise)
tpc_data$performance <- tpc_data$performance + rnorm(nrow(tpc_data), 0, 0.05)
tpc_data$performance[tpc_data$performance < 0] <- 0.001

# Visualize raw data
ggplot(tpc_data, aes(temp, performance)) +
 geom_point(size = 3) +
 labs(title = "Raw TPC Data",
      x = "Temperature (°C)",
      y = "Performance (development rate)") +
 theme_minimal()

# === STEP 3: Fit TPC models using rTPC ===

# Try multiple model types to find best fit
# Common models: Sharpe-Schoolfield, Briere, Gaussian, etc.

# Model 1: Sharpe-Schoolfield (high temperature inactivation)
fit_sharpeschoolhigh <- nls_multstart(
 performance ~ sharpeschoolhigh_1981(temp = temp, r_tref, e, eh, th, tref = 25),
 data = tpc_data,
 iter = 500,
 start_lower = c(r_tref = 0.1, e = 0.3, eh = 0.5, th = 285),
 start_upper = c(r_tref = 2, e = 1, eh = 3, th = 320),
 supp_errors = 'Y',
 convergence_count = 100
)

# Model 2: Briere (simpler, often good for insects)
fit_briere <- nls_multstart(
 performance ~ briere2_1999(temp = temp, tmin, tmax, a, b),
 data = tpc_data,
 iter = 500,
 start_lower = c(tmin = 5, tmax = 35, a = 0.0001, b = 1),
 start_upper = c(tmin = 15, tmax = 45, a = 0.01, b = 3),
 supp_errors = 'Y',
 convergence_count = 100
)

# Model 3: Gaussian (simple bell curve)
fit_gaussian <- nls_multstart(
 performance ~ gaussian_1987(temp = temp, rmax, topt, a),
 data = tpc_data,
 iter = 500,
 start_lower = c(rmax = 0.5, topt = 20, a = 5),
 start_upper = c(rmax = 1.5, topt = 30, a = 20),
 supp_errors = 'Y',
 convergence_count = 100
)

# === STEP 4: Compare model fits ===

# Calculate AIC for model selection
aic_sharp <- AIC(fit_sharpeschoolhigh)
aic_briere <- AIC(fit_briere)
aic_gaussian <- AIC(fit_gaussian)

cat("\nModel comparison (lower AIC is better):\n")
cat("Sharpe-Schoolfield AIC:", round(aic_sharp, 2), "\n")
cat("Briere AIC:", round(aic_briere, 2), "\n")
cat("Gaussian AIC:", round(aic_gaussian, 2), "\n")

# Select best model
best_model <- which.min(c(aic_sharp, aic_briere, aic_gaussian))
model_names <- c("Sharpe-Schoolfield", "Briere", "Gaussian")
cat("\nBest model:", model_names[best_model], "\n")

# === STEP 5: Visualize fitted models ===

# Generate predictions
temp_seq <- seq(5, 45, length.out = 200)

pred_sharp <- data.frame(
 temp = temp_seq,
 performance = predict(fit_sharpeschoolhigh, 
                       newdata = data.frame(temp = temp_seq)),
 model = "Sharpe-Schoolfield"
)

pred_briere <- data.frame(
 temp = temp_seq,
 performance = predict(fit_briere, 
                       newdata = data.frame(temp = temp_seq)),
 model = "Briere"
)

pred_gaussian <- data.frame(
 temp = temp_seq,
 performance = predict(fit_gaussian, 
                       newdata = data.frame(temp = temp_seq)),
 model = "Gaussian"
)

all_preds <- bind_rows(pred_sharp, pred_briere, pred_gaussian)

ggplot() +
 geom_line(data = all_preds, aes(temp, performance, color = model), size = 1) +
 geom_point(data = tpc_data, aes(temp, performance), size = 3) +
 labs(title = "Fitted TPC Models",
      x = "Temperature (°C)",
      y = "Performance") +
 theme_minimal() +
 theme(legend.position = "bottom")

# === STEP 6: Extract parameters from best model ===
# Using Sharpe-Schoolfield as example

params <- coef(fit_sharpeschoolhigh)
print(params)

# Get thermal optimum and limits
thermal_traits <- get_thermallimits(fit_sharpeschoolhigh)
cat("\nThermal traits:\n")
print(thermal_traits)

# === STEP 7: Create spatial temperature data ===
# Option A: Load real climate data (e.g., WorldClim)
# library(geodata)
# temp_raster <- worldclim_global(var = "tavg", res = 10, path = tempdir())

# Option B: Create example temperature raster
temp_raster <- rast(ncols = 200, nrows = 200, 
                    xmin = -10, xmax = 10, 
                    ymin = -10, ymax = 10,
                    crs = "EPSG:4326")

# Simulate realistic temperature gradient
values(temp_raster) <- 15 + 10 * (1:ncell(temp_raster) / ncell(temp_raster)) + 
 rnorm(ncell(temp_raster), 0, 2)

names(temp_raster) <- "temperature"

plot(temp_raster, main = "Mean Annual Temperature (°C)",
     col = hcl.colors(100, "RdYlBu", rev = TRUE))

# === STEP 8: Convert TPC to Probability of Presence ===

tpc_to_probability <- function(temp_raster, fitted_model, 
                               model_type = "sharpeschoolhigh",
                               normalize = TRUE) {
 
 # Extract temperature values
 temp_vals <- values(temp_raster, mat = FALSE)
 
 # Create prediction data frame
 pred_data <- data.frame(temp = temp_vals)
 
 # Predict performance
 performance <- predict(fitted_model, newdata = pred_data)
 
 # Handle NAs and negative values
 performance[is.na(performance)] <- 0
 performance[performance < 0] <- 0
 
 # Normalize to 0-1 if requested
 if (normalize) {
  max_perf <- max(performance, na.rm = TRUE)
  if (max_perf > 0) {
   probability <- performance / max_perf
  } else {
   probability <- performance
  }
 } else {
  probability <- performance
 }
 
 # Create output raster
 prob_raster <- temp_raster
 values(prob_raster) <- probability
 names(prob_raster) <- "prob_presence"
 
 return(prob_raster)
}

# Apply conversion using best fitted model
prob_map <- tpc_to_probability(temp_raster, 
                               fitted_model = fit_sharpeschoolhigh,
                               normalize = TRUE)

# === STEP 9: Visualize probability map ===

plot(prob_map, 
     main = "Probability of Presence\n(based on thermal performance)",
     col = hcl.colors(100, "YlOrRd", rev = TRUE),
     range = c(0, 1))

# === STEP 10: Create multiple probability maps with uncertainty ===

# Bootstrap confidence intervals (advanced)
bootstrap_tpc_predictions <- function(temp_raster, tpc_data, model_formula,
                                      n_bootstrap = 100) {
 
 temp_vals <- values(temp_raster, mat = FALSE)
 n_cells <- length(temp_vals)
 prob_matrix <- matrix(NA, nrow = n_cells, ncol = n_bootstrap)
 
 for (i in 1:n_bootstrap) {
  # Resample data
  boot_indices <- sample(1:nrow(tpc_data), replace = TRUE)
  boot_data <- tpc_data[boot_indices, ]
  
  # Fit model
  tryCatch({
   boot_fit <- nls_multstart(
    performance ~ sharpeschoolhigh_1981(temp = temp, r_tref, e, eh, th, tref = 25),
    data = boot_data,
    iter = 100,
    start_lower = c(r_tref = 0.1, e = 0.3, eh = 0.5, th = 285),
    start_upper = c(r_tref = 2, e = 1, eh = 3, th = 320),
    supp_errors = 'Y',
    convergence_count = 50
   )
   
   # Predict
   boot_perf <- predict(boot_fit, newdata = data.frame(temp = temp_vals))
   boot_perf[boot_perf < 0] <- 0
   boot_perf[is.na(boot_perf)] <- 0
   prob_matrix[, i] <- boot_perf / max(boot_perf, na.rm = TRUE)
  }, error = function(e) {
   # If fit fails, skip this iteration
   return(NULL)
  })
  
  if (i %% 10 == 0) cat("Bootstrap iteration", i, "of", n_bootstrap, "\n")
 }
 
 # Calculate mean and CI
 mean_prob <- rowMeans(prob_matrix, na.rm = TRUE)
 lower_ci <- apply(prob_matrix, 1, quantile, probs = 0.025, na.rm = TRUE)
 upper_ci <- apply(prob_matrix, 1, quantile, probs = 0.975, na.rm = TRUE)
 uncertainty <- upper_ci - lower_ci
 
 # Create raster stack
 prob_mean <- temp_raster
 values(prob_mean) <- mean_prob
 names(prob_mean) <- "mean_probability"
 
 prob_lower <- temp_raster
 values(prob_lower) <- lower_ci
 names(prob_lower) <- "lower_ci"
 
 prob_upper <- temp_raster
 values(prob_upper) <- upper_ci
 names(prob_upper) <- "upper_ci"
 
 prob_uncertainty <- temp_raster
 values(prob_uncertainty) <- uncertainty
 names(prob_uncertainty) <- "uncertainty"
 
 result <- c(prob_mean, prob_lower, prob_upper, prob_uncertainty)
 
 return(result)
}

# Run bootstrap (this takes a few minutes)
# Uncomment to run:
# prob_stack_bootstrap <- bootstrap_tpc_predictions(temp_raster, tpc_data, n_bootstrap = 50)
# plot(prob_stack_bootstrap)

# === STEP 11: Threshold for binary presence/absence ===

# Method 1: Fixed threshold
threshold <- 0.5
binary_map <- prob_map >= threshold

plot(binary_map, 
     main = paste0("Predicted Presence (threshold = ", threshold, ")"),
     col = c("white", "darkgreen"),
     legend = FALSE)

# Method 2: Percentile-based threshold (e.g., top 20% suitable)
threshold_percentile <- quantile(values(prob_map), probs = 0.80, na.rm = TRUE)
binary_map_percentile <- prob_map >= threshold_percentile

# === STEP 12: Calculate habitat suitability statistics ===

# Total suitable area
suitable_cells <- sum(values(binary_map) == 1, na.rm = TRUE)
total_cells <- sum(!is.na(values(binary_map)))
percent_suitable <- (suitable_cells / total_cells) * 100

cat("\n=== Habitat Suitability Summary ===\n")
cat("Threshold:", threshold, "\n")
cat("Suitable cells:", suitable_cells, "\n")
cat("Percent of area suitable:", round(percent_suitable, 2), "%\n")

# Probability distribution
prob_vals <- values(prob_map, mat = FALSE)
cat("\nProbability statistics:\n")
cat("Mean:", round(mean(prob_vals, na.rm = TRUE), 3), "\n")
cat("Median:", round(median(prob_vals, na.rm = TRUE), 3), "\n")
cat("Range:", round(min(prob_vals, na.rm = TRUE), 3), "-", 
    round(max(prob_vals, na.rm = TRUE), 3), "\n")

# Histogram of probabilities
hist(prob_vals, breaks = 50, 
     main = "Distribution of Presence Probabilities",
     xlab = "Probability of Presence",
     col = "steelblue")

# === STEP 13: Save outputs ===

# Save rasters
writeRaster(prob_map, "probability_presence.tif", overwrite = TRUE)
writeRaster(binary_map, "binary_presence.tif", overwrite = TRUE)

# Save model parameters
model_summary <- data.frame(
 parameter = names(params),
 value = as.numeric(params),
 model = "Sharpe-Schoolfield"
)

write.csv(model_summary, "tpc_parameters.csv", row.names = FALSE)

# Save thermal traits
write.csv(as.data.frame(thermal_traits), "thermal_traits.csv", row.names = FALSE)

cat("\n=== Analysis complete! ===\n")
cat("Files saved:\n")
cat("- probability_presence.tif\n")
cat("- binary_presence.tif\n")
cat("- tpc_parameters.csv\n")
cat("- thermal_traits.csv\n")

# === BONUS: Side-by-side comparison plot ===

par(mfrow = c(2, 2))
plot(temp_raster, main = "Temperature (°C)", 
     col = hcl.colors(100, "RdYlBu", rev = TRUE))
plot(prob_map, main = "Probability of Presence",
     col = hcl.colors(100, "YlOrRd", rev = TRUE))
plot(binary_map, main = "Binary Presence",
     col = c("white", "darkgreen"))
hist(values(prob_map), breaks = 30, main = "Probability Distribution",
     xlab = "Probability", col = "steelblue")
