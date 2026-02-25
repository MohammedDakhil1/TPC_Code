
#####################Complete Workflow: From Trait Data → TPC → Probability Map
# === PACKAGES ===
library(rTPC)
library(nls.multstart)
library(terra)
library(ggplot2)
library(tidyverse)

# ==========================================
# STEP 1: Prepare Your Trait-Temperature Data
# ==========================================

# Example: You measured body size, growth rate, or survival at different temps
# Your data structure:
trait_temp_data <- data.frame(
 temperature = c(10, 10, 10, 12, 12, 12, 15, 15, 15, 18, 18, 18, 20, 20, 20,
                 22, 22, 22, 25, 25, 25, 27, 27, 27, 30, 30, 30, 32, 32, 32,
                 35, 35, 35, 37, 37, 37, 40, 40, 40),
 trait_value = c(2.1, 2.3, 2.0, 3.5, 3.2, 3.4, 5.1, 5.4, 5.0, 7.2, 6.8, 7.0, 8.5, 8.9, 8.7,
                 9.2, 9.5, 9.4, 9.8, 9.6, 9.7, 9.3, 9.1, 9.0, 7.8, 8.2, 8.0, 6.1, 5.8, 6.0,
                 3.5, 3.8, 3.6, 1.2, 1.5, 1.3, 0.2, 0.3, 0.1)
)
# If loading from file:
# trait_temp_data <- read.csv("my_trait_data.csv")

# Visualize raw data
ggplot(trait_temp_data, aes(temperature, trait_value)) +
 geom_point(size = 3, alpha = 0.6) +
 labs(title = "Raw Trait-Temperature Data",
      x = "Temperature (°C)",
      y = "Trait Value (e.g., growth rate, body size)") +
 theme_minimal()

# ==========================================
# STEP 2: Calculate Mean Trait per Temperature
# ==========================================
# (if you have replicates)
trait_summary <- trait_temp_data %>%
 group_by(temperature) %>%
 summarise(
  mean_trait = mean(trait_value, na.rm = TRUE),
  se_trait = sd(trait_value, na.rm = TRUE) / sqrt(n()),
  n = n()
 )

print(trait_summary)

# Visualize with error bars
ggplot(trait_summary, aes(temperature, mean_trait)) +
 geom_point(size = 3) +
 geom_errorbar(aes(ymin = mean_trait - se_trait, 
                   ymax = mean_trait + se_trait), 
               width = 0.5) +
 labs(title = "Mean Trait Values Across Temperatures",
      x = "Temperature (°C)",
      y = "Mean Trait Value") +
 theme_minimal()

# ==========================================
# STEP 3: Standardize/Transform Trait if Needed
# ==========================================

# Option A: Normalize to 0-1 scale (relative performance)
trait_summary <- trait_summary %>%
 mutate(performance = mean_trait / max(mean_trait))
trait_summary
# Option B: Convert to rate (if measuring time to development)
# If trait is "days to maturity", convert to rate:
# trait_summary$performance <- 1 / trait_summary$mean_trait

# Option C: Log-transform if needed
# trait_summary$performance <- log(trait_summary$mean_trait)

# For this example, use normalized values
ggplot(trait_summary, aes(temperature, performance)) +
 geom_point(size = 3) +
 labs(title = "Normalized Performance (TPC)",
      x = "Temperature (°C)",
      y = "Relative Performance (0-1)") +
 theme_minimal()

# ==========================================
# STEP 4: Fit TPC Models
# ==========================================

# Prepare data for modeling
tpc_fit_data <- data.frame(
 temp = trait_summary$temperature,
 performance = trait_summary$performance
)

# Remove any NA values
tpc_fit_data <- na.omit(tpc_fit_data)
tpc_fit_data
# === Try Multiple TPC Models ===

# MODEL 1: Sharpe-Schoolfield (most common, accounts for high-temp decline)
cat("\n=== Fitting Sharpe-Schoolfield Model ===\n")

fit_sharpe <- nls_multstart(
 performance ~ sharpeschoolhigh_1981(temp = temp, r_tref, e, eh, th, tref = 20),
 data = tpc_fit_data,
 iter = 500,
 start_lower = c(r_tref = 0.01, e = 0.1, eh = 0.2, th = 273 + 25),
 start_upper = c(r_tref = 5, e = 2, eh = 10, th = 273 + 45),
 supp_errors = 'Y',
 convergence_count = 100
)

if (!is.null(fit_sharpe)) {
 cat("✓ Sharpe-Schoolfield converged\n")
 print(summary(fit_sharpe))
}

# MODEL 2: Briere (good for insects, asymmetric curve)
cat("\n=== Fitting Briere Model ===\n")

fit_briere <- nls_multstart(
 performance ~ briere2_1999(temp = temp, tmin, tmax, a, b),
 data = tpc_fit_data,
 iter = 500,
 start_lower = c(tmin = 0, tmax = 30, a = 0.00001, b = 0.1),
 start_upper = c(tmin = 20, tmax = 50, a = 0.1, b = 5),
 supp_errors = 'Y',
 convergence_count = 100
)

if (!is.null(fit_briere)) {
 cat("✓ Briere converged\n")
 print(summary(fit_briere))
}

# MODEL 3: Gaussian (simple bell curve)
cat("\n=== Fitting Gaussian Model ===\n")

fit_gaussian <- nls_multstart(
 performance ~ gaussian_1987(temp = temp, rmax, topt, a),
 data = tpc_fit_data,
 iter = 500,
 start_lower = c(rmax = 0.1, topt = 15, a = 1),
 start_upper = c(rmax = 2, topt = 35, a = 30),
 supp_errors = 'Y',
 convergence_count = 100
)

if (!is.null(fit_gaussian)) {
 cat("✓ Gaussian converged\n")
 print(summary(fit_gaussian))
}

# MODEL 4: Quadratic (simplest, often good starting point)
cat("\n=== Fitting Quadratic Model ===\n")

fit_quadratic <- nls_multstart(
 performance ~ quadratic_2008(temp = temp, a, b, c),
 data = tpc_fit_data,
 iter = 500,
 start_lower = c(a = -1, b = -1, c = -10),
 start_upper = c(a = 1, b = 5, c = 50),
 supp_errors = 'Y',
 convergence_count = 100
)

if (!is.null(fit_quadratic)) {
 cat("✓ Quadratic converged\n")
 print(summary(fit_quadratic))
}

# ==========================================
# STEP 5: Compare Models Using AIC
# ==========================================

# Collect all successful fits
fitted_models <- list()
model_names <- c()

if (!is.null(fit_sharpe)) {
 fitted_models$sharpe <- fit_sharpe
 model_names <- c(model_names, "Sharpe-Schoolfield")
}
if (!is.null(fit_briere)) {
 fitted_models$briere <- fit_briere
 model_names <- c(model_names, "Briere")
}
if (!is.null(fit_gaussian)) {
 fitted_models$gaussian <- fit_gaussian
 model_names <- c(model_names, "Gaussian")
}
if (!is.null(fit_quadratic)) {
 fitted_models$quadratic <- fit_quadratic
 model_names <- c(model_names, "Quadratic")
}

# Calculate AIC for each model
aic_values <- sapply(fitted_models, AIC)
r2_values <- sapply(fitted_models, function(m) {
 1 - sum(residuals(m)^2) / sum((tpc_fit_data$performance - mean(tpc_fit_data$performance))^2)
})

model_comparison <- data.frame(
 Model = model_names,
 AIC = round(aic_values, 3),
 R_squared = round(r2_values, 3),
 Delta_AIC = round(aic_values - min(aic_values), 3)
)

model_comparison <- model_comparison %>%
 arrange(AIC)

cat("\n=== MODEL COMPARISON ===\n")
print(model_comparison)

# Select best model
best_model_name <- model_comparison$Model[1]
best_model <- fitted_models[[which(model_names == best_model_name)]]

cat("\n✓ Best model:", best_model_name, 
    "(lowest AIC =", round(min(aic_values), 2), ")\n")

# ==========================================
# STEP 6: Visualize All Fitted Models
# ==========================================

# Generate smooth predictions
temp_predict <- seq(min(tpc_fit_data$temp) - 5, 
                    max(tpc_fit_data$temp) + 5, 
                    length.out = 200)

# Predict for all models
all_predictions <- data.frame(temp = temp_predict)

for (i in seq_along(fitted_models)) {
 model <- fitted_models[[i]]
 model_name <- names(fitted_models)[i]
 
 tryCatch({
  pred <- predict(model, newdata = data.frame(temp = temp_predict))
  all_predictions[[model_name]] <- pred
 }, error = function(e) {
  all_predictions[[model_name]] <- NA
 })
}

# Reshape for plotting
predictions_long <- all_predictions %>%
 pivot_longer(-temp, names_to = "model", values_to = "performance")

# Plot
ggplot() +
 geom_line(data = predictions_long, 
           aes(temp, performance, color = model), 
           size = 1, alpha = 0.7) +
 geom_point(data = tpc_fit_data, 
            aes(temp, performance), 
            size = 3, color = "black") +
 labs(title = "TPC Model Comparison",
      subtitle = paste("Best model:", best_model_name),
      x = "Temperature (°C)",
      y = "Performance",
      color = "Model") +
 theme_minimal() +
 theme(legend.position = "bottom")

#ggsave("tpc_model_comparison.png", width = 10, height = 6, dpi = 300)

# ==========================================
# STEP 7: Extract Thermal Parameters from Best Model, see next code solving the error
# ==========================================

# Get thermal limits and optimum
#thermal_params <- get_thermallimits(best_model)

#cat("\n=== THERMAL PARAMETERS (from best model) ===\n")
#print(as.data.frame(thermal_params))

# Additional metrics
params <- coef(best_model)
cat("\nModel coefficients:\n")
print(params)

# Calculate optimum manually for some models
temp_fine <- seq(0, 50, by = 0.1)
perf_fine <- predict(best_model, newdata = data.frame(temp = temp_fine))
perf_fine[perf_fine < 0] <- 0

t_opt <- temp_fine[which.max(perf_fine)]
max_perf <- max(perf_fine, na.rm = TRUE)

cat("\nThermal optimum (Topt):", round(t_opt, 2), "°C\n")
cat("Maximum performance:", round(max_perf, 3), "\n")

# Critical thermal limits (where performance drops to 50% or 10%)
perf_50 <- max_perf * 0.5
perf_10 <- max_perf * 0.1

ct_min_50 <- temp_fine[which(perf_fine >= perf_50)[1]]
ct_max_50 <- temp_fine[tail(which(perf_fine >= perf_50), 1)]
ct_min_10 <- temp_fine[which(perf_fine >= perf_10)[1]]
ct_max_10 <- temp_fine[tail(which(perf_fine >= perf_10), 1)]

cat("\nCritical thermal limits (50% performance):", 
    round(ct_min_50, 2), "-", round(ct_max_50, 2), "°C\n")
cat("Critical thermal limits (10% performance):", 
    round(ct_min_10, 2), "-", round(ct_max_10, 2), "°C\n")

# Thermal breadth
thermal_breadth_50 <- ct_max_50 - ct_min_50
thermal_breadth_10 <- ct_max_10 - ct_min_10

cat("Thermal breadth (50%):", round(thermal_breadth_50, 2), "°C\n")
cat("Thermal breadth (10%):", round(thermal_breadth_10, 2), "°C\n")

# ==========================================
# STEP 8: Create Spatial Probability Map
# ==========================================

# Load or create temperature raster
# Option 1: Load real data
# library(geodata)
# temp_raster <- worldclim_global(var = "tavg", res = 10, path = tempdir())
# temp_raster <- mean(temp_raster) # if you have 12 monthly layers

# Option 2: Create example data
temp_raster <- rast(ncols = 300, nrows = 300,
                    xmin = -120, xmax = -70,
                    ymin = 25, ymax = 50,
                    crs = "EPSG:4326")

# Simulate temperature gradient (cooler north, warmer south)
coords <- crds(temp_raster)
values(temp_raster) <- 5 + (50 - coords[,2]) * 0.5 + rnorm(nrow(coords), 0, 2)
names(temp_raster) <- "temperature"

plot(temp_raster, main = "Temperature (°C)",
     col = hcl.colors(100, "RdYlBu", rev = TRUE))

# Convert TPC to probability
tpc_to_prob_map <- function(temp_raster, fitted_model, normalize = TRUE) {
 
 temp_vals <- values(temp_raster, mat = FALSE)
 temp_vals
 # Predict performance
 pred_data <- data.frame(temp = temp_vals)
 performance <- predict(fitted_model, newdata = pred_data)
 
 # Clean up predictions
 performance[is.na(performance)] <- 0
 performance[performance < 0] <- 0
 performance[is.infinite(performance)] <- 0
 
 # Normalize to 0-1
 if (normalize) {
  max_perf <- max(performance, na.rm = TRUE)
  if (max_perf > 0) {
   probability <- performance / max_perf
  } else {
   probability <- rep(0, length(performance))
  }
 } else {
  probability <- performance
 }
 
 # Create raster
 prob_raster <- temp_raster
 values(prob_raster) <- probability
 names(prob_raster) <- "probability_presence"
 
 return(prob_raster)
}

# Generate probability map
prob_map <- tpc_to_prob_map(temp_raster, best_model, normalize = TRUE)

plot(prob_map, 
     main = paste("Probability of Presence\n(Based on", best_model_name, "TPC)"),
     col = hcl.colors(100, "YlOrRd", rev = TRUE),
     range = c(0, 1))
####################################Modified code using NicheMapR Ectotherm function (ectotherm or microclimate model)
# ==========================================
# STEP 8: Create Spatial Probability Map using NicheMapR
# ==========================================

#################################Fixed error, modified2
folder <- climate_folder
save(folder, file = paste0(.libPaths()[1], "/gcfolder.rda"))

gcfolder <- paste(.libPaths()[1],"/gcfolder.rda",sep="")
if(file.exists(gcfolder) == FALSE){
 folder <- "c:/globalclimate"
 if(file.exists(paste0(folder,"/global_climate.nc")) == FALSE){
  message("You don't appear to have the global climate data set...")
  stop()
 }
}else{
 load(gcfolder)
}
# ==========================================
# ==========================================
# STEP 8: Create Spatial Probability Map using NicheMapR, take long time extracting global temp data
# ==========================================

library(NicheMapR)
library(terra)

# ------------------------------------------
# 8a: Download & register the global climate dataset (ONE-TIME SETUP)
# ------------------------------------------
# micro_global() looks for climate data in two ways:
#   1. It checks .libPaths()[1]/gcfolder.rda — a small file that stores
#      the PATH to your climate data folder
#   2. If that file doesn't exist, it falls back to "c:/globalclimate"
#
# The error occurs because NEITHER of these is working on your system.
# Below are three approaches — pick ONE that works for you.

# --- APPROACH 1 (Recommended): Download data + manually create the pointer ---
# Choose where to store the ~0.5 GB climate data
climate_folder <- "~/NicheMapR_climate"  # Change to any folder you want

# Create the folder
if (!dir.exists(climate_folder)) {
 dir.create(climate_folder, recursive = TRUE)
}

# Download the data (only needs to run ONCE — comment out after success)
# When prompted, type 'y' and press Enter
get.global.climate(folder = climate_folder)

# CRITICAL FIX: Manually create the pointer file that micro_global looks for.
# get.global.climate() is supposed to create this, but it sometimes fails
# due to write permissions in the R library folder.
# This saves a small .rda file that tells micro_global() where your data lives.

folder <- climate_folder  # this variable name must be exactly "folder"
save(folder, file = paste0(.libPaths()[1], "/gcfolder.rda"))

# Verify the pointer was saved:
cat("Pointer file location:", paste0(.libPaths()[1], "/gcfolder.rda"), "\n")
cat("Pointer file exists:", file.exists(paste0(.libPaths()[1], "/gcfolder.rda")), "\n")

# Verify the climate data files exist:
cat("global_climate.nc exists:", 
    file.exists(paste0(climate_folder, "/global_climate.nc")), "\n")
cat("MODIS.nc exists:", 
    file.exists(paste0(climate_folder, "/MODIS.nc")), "\n")


# --- APPROACH 2: If you can't write to .libPaths(), use the default path ---
# micro_global() hardcodes "c:/globalclimate" as a fallback (Windows).
# On Mac/Linux, create a symlink to trick it:
#
# On Windows:
#   get.global.climate(folder = "c:/globalclimate")
#
# On Mac/Linux (create a symlink to the default Windows path):
#   get.global.climate(folder = "~/NicheMapR_climate")
#   Then in R:
#   file.symlink("~/NicheMapR_climate", "/tmp/globalclimate")
#   And override as in Approach 1.


# --- APPROACH 3: Skip micro_global entirely, use micro_ncep instead ---
# micro_ncep() does NOT need the pre-downloaded global climate dataset.
# It fetches NCEP reanalysis data directly from the internet for specific
# date ranges. This is actually better for many use cases because it gives
# you real historical weather rather than long-term averages.
# See section 8c_alternative below.


# ------------------------------------------
# 8b: Define your study region as a grid of coordinates
# ------------------------------------------
lon_seq <- seq(-120, -70, by = 2)
lat_seq <- seq(25, 50, by = 2)
grid_points <- expand.grid(lon = lon_seq, lat = lat_seq)

cat("Total grid points to model:", nrow(grid_points), "\n")


# ------------------------------------------
# 8c: Run micro_global at each grid point (using Approach 1 or 2)
# ------------------------------------------
# Now that the pointer file exists, micro_global() should find the data.

mean_temps <- rep(NA_real_, nrow(grid_points))

for (i in seq_len(nrow(grid_points))) {
 
 cat("\rProcessing point", i, "of", nrow(grid_points), "...")
 
 tryCatch({
  micro <- micro_global(
   loc          = c(grid_points$lon[i], grid_points$lat[i]),
   Usrhyt       = 0.01,
   minshade     = 0,
   maxshade     = 90,
   timeinterval = 12,
   run.gads     = 2
  )
  
  metout <- as.data.frame(micro$metout)
  mean_temps[i] <- mean(metout$TALOC, na.rm = TRUE)
  
 }, error = function(e) {
  mean_temps[i] <<- NA
  cat(" [failed:", conditionMessage(e), "]")
 })
}

cat("\nDone. Successfully modeled", sum(!is.na(mean_temps)), 
    "of", nrow(grid_points), "locations.\n")


# ------------------------------------------
# 8c_alternative: Use micro_ncep instead (Approach 3 — no download needed)
# ------------------------------------------
# micro_ncep fetches NCEP reanalysis weather data from the internet
# for a specific date range. No local climate database required.
# Slower per point (web requests), but no setup headaches.

# Uncomment below to use this approach instead of 8c:

# mean_temps <- rep(NA_real_, nrow(grid_points))
# 
# for (i in seq_len(nrow(grid_points))) {
#   
#   cat("\rProcessing point", i, "of", nrow(grid_points), "...")
#   
#   tryCatch({
#     # Fetch real NCEP weather for 2023
#     micro <- micro_ncep(
#       loc      = c(grid_points$lon[i], grid_points$lat[i]),
#       dstart   = "01/01/2023",
#       dfinish  = "31/12/2023",
#       Usrhyt   = 0.01,
#       minshade = 0,
#       maxshade = 90
#     )
#     
#     metout <- as.data.frame(micro$metout)
#     mean_temps[i] <- mean(metout$TALOC, na.rm = TRUE)
#     
#   }, error = function(e) {
#     mean_temps[i] <<- NA
#     cat(" [failed:", conditionMessage(e), "]")
#   })
# }
# 
# cat("\nDone. Successfully modeled", sum(!is.na(mean_temps)),
#     "of", nrow(grid_points), "locations.\n")


# ------------------------------------------
# 8d: Convert point temperatures into a raster
# ------------------------------------------
grid_points$mean_temp <- mean_temps

temp_raster <- rast(
 xmin = min(lon_seq), xmax = max(lon_seq),
 ymin = min(lat_seq), ymax = max(lat_seq),
 res  = 2,
 crs  = "EPSG:4326"
)

temp_raster <- rasterize(
 x     = vect(grid_points, geom = c("lon", "lat"), crs = "EPSG:4326"),
 y     = temp_raster,
 field = "mean_temp"
)

names(temp_raster) <- "temperature"

plot(temp_raster,
     main = "NicheMapR Microclimate Temperature (°C)\n(mean annual at 1 cm height)",
     col  = hcl.colors(100, "RdYlBu", rev = TRUE))


# ------------------------------------------
# 8e: Apply the TPC to create a probability map
# ------------------------------------------
tpc_to_prob_map <- function(temp_raster, fitted_model, normalize = TRUE) {
 
 temp_vals <- values(temp_raster, mat = FALSE)
 performance <- predict(fitted_model, newdata = data.frame(temp = temp_vals))
 
 performance[is.na(performance)]       <- 0
 performance[performance < 0]          <- 0
 performance[is.infinite(performance)] <- 0
 
 if (normalize && max(performance, na.rm = TRUE) > 0) {
  probability <- performance / max(performance, na.rm = TRUE)
 } else {
  probability <- performance
 }
 
 prob_raster <- temp_raster
 values(prob_raster) <- probability
 names(prob_raster) <- "probability_presence"
 
 return(prob_raster)
}

prob_map <- tpc_to_prob_map(temp_raster, best_model, normalize = TRUE)
prob_map
plot(prob_map,
     main  = paste("Probability of Presence\n(NicheMapR +", best_model_name, "TPC)"),
     col   = hcl.colors(100, "YlOrRd", rev = TRUE),
     range = c(0, 1))
# ==========================================
# STEP 9: Create Binary Presence Map
# ==========================================

# Multiple threshold approaches
threshold_50 <- 0.5  # 50% of optimal performance
threshold_percentile <- quantile(values(prob_map), 0.75, na.rm = TRUE)

binary_map_50 <- prob_map >= threshold_50
binary_map_75th <- prob_map >= threshold_percentile

par(mfrow = c(1, 2))
plot(binary_map_50, main = "Presence (threshold = 0.5)",
     col = c("white", "darkgreen"), legend = FALSE)
plot(binary_map_75th, main = "Presence (75th percentile)",
     col = c("white", "darkgreen"), legend = FALSE)

# ==========================================
# STEP 10: Calculate Suitability Statistics
# ==========================================

suitable_area_50 <- sum(values(binary_map_50) == 1, na.rm = TRUE)
total_area <- sum(!is.na(values(binary_map_50)))
percent_suitable_50 <- (suitable_area_50 / total_area) * 100

cat("\n=== HABITAT SUITABILITY SUMMARY ===\n")
cat("Total pixels:", total_area, "\n")
cat("Suitable pixels (p > 0.5):", suitable_area_50, "\n")
cat("Percent suitable:", round(percent_suitable_50, 2), "%\n")

# Probability distribution
prob_values <- values(prob_map, mat = FALSE)
cat("\nProbability distribution:\n")
cat("  Mean:", round(mean(prob_values, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(prob_values, na.rm = TRUE), 3), "\n")
cat("  SD:", round(sd(prob_values, na.rm = TRUE), 3), "\n")
cat("  Range:", round(min(prob_values, na.rm = TRUE), 3), "-",
    round(max(prob_values, na.rm = TRUE), 3), "\n")

# ==========================================
# STEP 11: Comprehensive Visualization
# ==========================================

par(mfrow = c(2, 3))

# 1. Original data
plot(tpc_fit_data$temp, tpc_fit_data$performance, 
     pch = 19, cex = 2,
     main = "Original Data", xlab = "Temperature (°C)", 
     ylab = "Performance")
lines(temp_predict, predict(best_model, data.frame(temp = temp_predict)),
      col = "red", lwd = 2)

# 2. Temperature map
plot(temp_raster, main = "Temperature", 
     col = hcl.colors(100, "RdYlBu", rev = TRUE))

# 3. Probability map
plot(prob_map, main = "Probability of Presence",
     col = hcl.colors(100, "YlOrRd", rev = TRUE))

# 4. Binary presence
plot(binary_map_50, main = "Binary Presence (p > 0.5)",
     col = c("white", "darkgreen"), legend = FALSE)

# 5. Probability histogram
hist(prob_values, breaks = 50, col = "steelblue",
     main = "Probability Distribution", xlab = "Probability")

# 6. Temperature vs Probability scatter
temp_vals <- values(temp_raster, mat = FALSE)
plot(temp_vals, prob_values, pch = ".", col = rgb(0,0,0,0.1),
     main = "Temp vs Probability", xlab = "Temperature (°C)",
     ylab = "Probability")


#####################

#################### STEP 7: Thermal Limit Function

# Enhanced thermal limits calculator with error handling
calculate_thermal_limits_robust <- function(fitted_model, 
                                            temp_range = c(-5, 55),
                                            thresholds = c(0.9, 0.8, 0.5, 0.1)) {
 
 # Generate temperature sequence
 temp_seq <- seq(temp_range[1], temp_range[2], by = 0.05)
 
 # Predict performance with error handling
 pred_perf <- tryCatch({
  predict(fitted_model, newdata = data.frame(temp = temp_seq))
 }, error = function(e) {
  warning("Prediction failed: ", e$message)
  return(rep(NA, length(temp_seq)))
 })
 
 # Clean predictions
 pred_perf[pred_perf < 0] <- 0
 pred_perf[is.na(pred_perf)] <- 0
 pred_perf[is.infinite(pred_perf)] <- 0
 
 if (all(pred_perf == 0)) {
  warning("All predictions are zero or invalid")
  return(NULL)
 }
 
 # Find thermal optimum
 topt_idx <- which.max(pred_perf)
 topt <- temp_seq[topt_idx]
 rmax <- pred_perf[topt_idx]
 
 # Initialize results list
 results_list <- list(
  topt = topt,
  rmax = rmax,
  temp_seq = temp_seq,
  performance = pred_perf
 )
 
 # Calculate limits for each threshold
 for (thresh in thresholds) {
  thresh_perf <- rmax * thresh
  above_thresh <- pred_perf >= thresh_perf
  
  if (any(above_thresh)) {
   temps_above <- temp_seq[above_thresh]
   ctmin <- min(temps_above)
   ctmax <- max(temps_above)
   breadth <- ctmax - ctmin
   
   thresh_label <- paste0("thresh_", round(thresh * 100))
   results_list[[paste0("ctmin_", thresh_label)]] <- ctmin
   results_list[[paste0("ctmax_", thresh_label)]] <- ctmax
   results_list[[paste0("breadth_", thresh_label)]] <- breadth
  } else {
   thresh_label <- paste0("thresh_", round(thresh * 100))
   results_list[[paste0("ctmin_", thresh_label)]] <- NA
   results_list[[paste0("ctmax_", thresh_label)]] <- NA
   results_list[[paste0("breadth_", thresh_label)]] <- NA
  }
 }
 
 # Calculate temperature at specific performance levels
 results_list$temp_at_50perf_lower <- temp_seq[which(pred_perf >= rmax * 0.5)[1]]
 results_list$temp_at_50perf_upper <- temp_seq[tail(which(pred_perf >= rmax * 0.5), 1)]
 
 # Calculate skewness (asymmetry of TPC)
 if (!is.na(results_list$ctmin_thresh_50) && !is.na(results_list$ctmax_thresh_50)) {
  skewness <- (topt - results_list$ctmin_thresh_50) / 
   (results_list$ctmax_thresh_50 - results_list$ctmin_thresh_50)
  results_list$skewness <- skewness
 } else {
  results_list$skewness <- NA
 }
 
 return(results_list)
}

# Use the robust function
thermal_limits <- calculate_thermal_limits_robust(best_model)

# Convert to data frame for easy viewing
thermal_df <- data.frame(
 Parameter = names(thermal_limits)[1:11],  # exclude temp_seq and performance
 Value = unlist(thermal_limits[1:11])
)

print(thermal_df)

# Plot with all thresholds
plot(thermal_limits$temp_seq, thermal_limits$performance, 
     type = "l", lwd = 2,
     main = "TPC with Multiple Thermal Thresholds",
     xlab = "Temperature (°C)", 
     ylab = "Performance")

# Add threshold lines
colors <- c("darkgreen", "green", "orange", "red")
thresh_levels <- c(0.9, 0.8, 0.5, 0.1)

for (i in seq_along(thresh_levels)) {
 abline(h = thermal_limits$rmax * thresh_levels[i], 
        col = colors[i], lty = 2)
 
 thresh_label <- paste0("thresh_", round(thresh_levels[i] * 100))
 ctmin <- thermal_limits[[paste0("ctmin_", thresh_label)]]
 ctmax <- thermal_limits[[paste0("ctmax_", thresh_label)]]
 
 if (!is.na(ctmin) && !is.na(ctmax)) {
  segments(ctmin, 0, ctmin, thermal_limits$rmax * thresh_levels[i], 
           col = colors[i], lty = 3)
  segments(ctmax, 0, ctmax, thermal_limits$rmax * thresh_levels[i], 
           col = colors[i], lty = 3)
 }
}

# Add data points
points(tpc_fit_data$temp, tpc_fit_data$performance, pch = 19, cex = 1.5)

# Add optimum
abline(v = thermal_limits$topt, col = "blue", lwd = 2, lty = 2)

legend("topright", 
       legend = c("90%", "80%", "50%", "10%", "Topt"),
       col = c(colors, "blue"),
       lty = 2, lwd = 1.5, cex = 0.7,
       title = "Performance threshold")

###########################Another long code for thermal limit
# ==========================================
# CORRECTED: Extract Thermal Parameters
# ==========================================

# Use the robust function
thermal_limits <- calculate_thermal_limits_robust(best_model)

# Convert to data frame - EXCLUDE temp_seq and performance vectors
thermal_scalar_names <- names(thermal_limits)[!names(thermal_limits) %in% c("temp_seq", "performance")]

thermal_df <- data.frame(
 Parameter = thermal_scalar_names,
 Value = unlist(thermal_limits[thermal_scalar_names])
)

cat("\n=== THERMAL PARAMETERS ===\n")
print(thermal_df)

# Extract key values for easy reference
topt <- thermal_limits$topt
rmax <- thermal_limits$rmax
ctmin_50 <- thermal_limits$ctmin_thresh_50
ctmax_50 <- thermal_limits$ctmax_thresh_50
breadth_50 <- thermal_limits$breadth_thresh_50

cat("\n=== KEY THERMAL TRAITS ===\n")
cat("Thermal optimum (Topt):", round(topt, 2), "°C\n")
cat("Maximum performance (rmax):", round(rmax, 3), "\n")
cat("Critical thermal minimum (CTmin, 50%):", round(ctmin_50, 2), "°C\n")
cat("Critical thermal maximum (CTmax, 50%):", round(ctmax_50, 2), "°C\n")
cat("Thermal breadth (50%):", round(breadth_50, 2), "°C\n")

# ==========================================
# Visualize Thermal Limits
# ==========================================

# Plot with all thresholds
plot(thermal_limits$temp_seq, thermal_limits$performance, 
     type = "l", lwd = 3, col = "black",
     main = "TPC with Multiple Thermal Thresholds",
     xlab = "Temperature (°C)", 
     ylab = "Performance",
     ylim = c(0, max(thermal_limits$performance, na.rm = TRUE) * 1.1))

# Add original data points
points(tpc_fit_data$temp, tpc_fit_data$performance, 
       pch = 19, cex = 2, col = "blue")

# Add threshold lines and limits
colors <- c("darkgreen", "green", "orange", "red")
thresh_levels <- c(0.9, 0.8, 0.5, 0.1)
thresh_names <- c("90%", "80%", "50%", "10%")

for (i in seq_along(thresh_levels)) {
 # Horizontal line at threshold
 abline(h = rmax * thresh_levels[i], 
        col = colors[i], lty = 2, lwd = 1.5)
 
 # Get thermal limits for this threshold
 thresh_label <- paste0("thresh_", round(thresh_levels[i] * 100))
 ctmin <- thermal_limits[[paste0("ctmin_", thresh_label)]]
 ctmax <- thermal_limits[[paste0("ctmax_", thresh_label)]]
 
 # Add vertical lines at limits
 if (!is.na(ctmin) && !is.na(ctmax)) {
  segments(ctmin, 0, ctmin, rmax * thresh_levels[i], 
           col = colors[i], lty = 3, lwd = 1.5)
  segments(ctmax, 0, ctmax, rmax * thresh_levels[i], 
           col = colors[i], lty = 3, lwd = 1.5)
  
  # Add shaded region for 50% threshold
  if (thresh_levels[i] == 0.5) {
   rect(ctmin, 0, ctmax, rmax * 0.5, 
        col = rgb(1, 0.5, 0, 0.15), border = NA)
  }
 }
}

# Add optimum line
abline(v = topt, col = "red", lwd = 2.5, lty = 2)
text(topt, max(thermal_limits$performance) * 1.05, 
     paste0("Topt = ", round(topt, 1), "°C"), 
     col = "red", cex = 0.9, font = 2)

# Add legend
legend("topright", 
       legend = c(thresh_names, "Topt", "Data"),
       col = c(colors, "red", "blue"),
       lty = c(rep(2, 4), 2, NA), 
       pch = c(rep(NA, 5), 19),
       lwd = c(rep(1.5, 4), 2.5, NA),
       cex = 0.8,
       title = "Performance threshold",
       bg = "white")

# ==========================================
# STEP 12: Save All Outputs
# ==========================================
setwd("C:/[[[[_Landscape_Ecology_Major_Revision/[[[[[[[[[[[[[[[[[_IOWA_US/[[[[[[[[[[[[[[[[[[[[[_University of Mary Washington/University of Illinois Urbana-Champaign/[[[[[[[[[[[[[[[[[_latest_CV/[[[[[[[[[[[_Notre_Dame_Biko_LAb/[[[[[[[[[[[_output")
# Save rasters
writeRaster(prob_map, "probability_map.tif", overwrite = TRUE)
writeRaster(binary_map_50, "binary_presence_50.tif", overwrite = TRUE)
writeRaster(temp_raster, "temperature_map.tif", overwrite = TRUE)

# Save TPC data
write.csv(tpc_fit_data, "tpc_fitted_data.csv", row.names = FALSE)

# Save model comparison
write.csv(model_comparison, "model_comparison.csv", row.names = FALSE)

# Save thermal parameters
thermal_summary <- data.frame(
 Parameter = c("T_optimum", "Max_Performance", "CT_min_50", "CT_max_50",
               "Thermal_Breadth_50", "CT_min_10", "CT_max_10", "Thermal_Breadth_10"),
 Value = c(t_opt, max_perf, ct_min_50, ct_max_50, thermal_breadth_50,
           ct_min_10, ct_max_10, thermal_breadth_10)
)
write.csv(thermal_summary, "thermal_parameters.csv", row.names = FALSE)

# Save best model coefficients
model_coefs <- data.frame(
 Model = best_model_name,
 Parameter = names(params),
 Value = as.numeric(params)
)
write.csv(model_coefs, "best_model_coefficients.csv", row.names = FALSE)

# Save final plot
png("complete_tpc_analysis.png", width = 12, height = 8, units = "in", res = 300)
par(mfrow = c(2, 3))
plot(tpc_fit_data$temp, tpc_fit_data$performance, pch = 19, cex = 2,
     main = "TPC Fit", xlab = "Temperature (°C)", ylab = "Performance")
lines(temp_predict, predict(best_model, data.frame(temp = temp_predict)),
      col = "red", lwd = 2)
plot(temp_raster, main = "Temperature", col = hcl.colors(100, "RdYlBu", rev = TRUE))
plot(prob_map, main = "Probability", col = hcl.colors(100, "YlOrRd", rev = TRUE))
plot(binary_map_50, main = "Presence", col = c("white", "darkgreen"), legend = FALSE)
hist(prob_values, breaks = 50, col = "steelblue", main = "Distribution", xlab = "Probability")
plot(temp_vals, prob_values, pch = ".", col = rgb(0,0,0,0.1),
     main = "Temp-Prob Relationship", xlab = "Temp", ylab = "Prob")
dev.off()

cat("\n✓ Analysis complete! All outputs saved.\n")


# ==========================================
