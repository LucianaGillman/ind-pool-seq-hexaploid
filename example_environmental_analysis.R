#!/usr/bin/env Rscript
# Example: Environmental Variable Transformation for Bromus auleticus
#
# This script demonstrates how to apply environmental variable transformation
# methodology to the Bromus auleticus dataset, following the approach used in
# studies like the Drosophila melanogaster environmental adaptation paper.

# Source the transformation functions
source("environmental_variable_transformation.R")

# Load required libraries
library(dplyr)

# Load metadata
cat("Loading metadata...\n")
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# Display basic information
cat("\nMetadata summary:\n")
print(summary(metadata))
cat("\nNumber of populations:", length(unique(metadata$pop)), "\n")
cat("Number of samples:", nrow(metadata), "\n")

# For this example, we'll simulate environmental data extraction
# In a real scenario, you would extract these from climate databases
# using the lat/lon coordinates (e.g., from WorldClim, CHELSA, etc.)

# Simulate environmental data for each unique population
cat("\n=== Simulating environmental data for populations ===\n")
cat("(In practice, extract from climate databases using lat/lon)\n\n")

pop_summary <- metadata %>%
  group_by(pop, lat, lon) %>%
  summarise(n_samples = n(), .groups = 'drop')

set.seed(42)
# Simulate environmental variables based on latitude
# (Southern latitudes generally have different climate patterns)
pop_environmental <- pop_summary %>%
  mutate(
    # Mean annual temperature (°C) - varies with latitude
    mean_annual_temp = 18 + (lat + 32) * 0.5 + rnorm(n(), 0, 1),
    
    # Annual rainfall (mm) - example values
    annual_rainfall = 1000 + (lat + 33) * 20 + rnorm(n(), 0, 100),
    
    # Mean wind speed (km/h) - example values
    mean_wind_speed = 15 + rnorm(n(), 0, 2)
  )

cat("Environmental data for populations:\n")
print(pop_environmental)

# Apply transformations using different methods

cat("\n\n=== METHOD 1: Median-split transformation ===\n")
cat("Variables above median = 1, below median = 0\n")
cat("This is commonly used to identify high vs. low environmental conditions\n\n")

result_median <- transform_environmental_variables(
  data = pop_environmental,
  variables = c("mean_annual_temp", "annual_rainfall", "mean_wind_speed"),
  method = "median_split",
  suffix = "_high_low"
)

cat("Transformed data (first 10 rows):\n")
print(result_median$data %>% select(pop, ends_with("_high_low")), n = 10)

cat("\nTransformation parameters:\n")
report_median <- generate_transformation_report(
  result_median$transformation_info,
  output_file = "transformation_report_median.csv"
)
print(report_median)

# Summary statistics for binary variables
cat("\n\nBinary variable distributions:\n")
for (var in c("mean_annual_temp_high_low", "annual_rainfall_high_low", "mean_wind_speed_high_low")) {
  cat(sprintf("\n%s:\n", var))
  cat(sprintf("  High (1): %d populations (%.1f%%)\n", 
              sum(result_median$data[[var]] == 1, na.rm = TRUE),
              100 * mean(result_median$data[[var]] == 1, na.rm = TRUE)))
  cat(sprintf("  Low (0): %d populations (%.1f%%)\n", 
              sum(result_median$data[[var]] == 0, na.rm = TRUE),
              100 * mean(result_median$data[[var]] == 0, na.rm = TRUE)))
}

cat("\n\n=== METHOD 2: Custom threshold transformation ===\n")
cat("Using biologically meaningful thresholds\n\n")

# Define custom thresholds based on biological relevance
custom_thresholds <- list(
  mean_annual_temp = 18,    # Temperature threshold (°C)
  annual_rainfall = 1000,   # Rainfall threshold (mm)
  mean_wind_speed = 15      # Wind speed threshold (km/h)
)

result_threshold <- transform_environmental_variables(
  data = pop_environmental,
  variables = c("mean_annual_temp", "annual_rainfall", "mean_wind_speed"),
  method = "threshold",
  thresholds = custom_thresholds,
  suffix = "_binary"
)

cat("Transformation parameters:\n")
report_threshold <- generate_transformation_report(
  result_threshold$transformation_info,
  output_file = "transformation_report_threshold.csv"
)
print(report_threshold)

# Save the transformed data
output_data <- result_threshold$data
write.csv(output_data, "population_environmental_binary.csv", row.names = FALSE)
cat("\n\nTransformed data saved to: population_environmental_binary.csv\n")

# Create a mapping between samples and their population's environmental classification
cat("\n\n=== Creating sample-level environmental classifications ===\n")

# Merge with original metadata
sample_environmental <- metadata %>%
  left_join(
    result_threshold$data %>% 
      select(pop, mean_annual_temp_binary, annual_rainfall_binary, mean_wind_speed_binary),
    by = "pop"
  )

cat("Sample data with environmental classifications (first 20 rows):\n")
print(head(sample_environmental, 20))

# Save sample-level data
write.csv(sample_environmental, "samples_with_environmental_binary.csv", row.names = FALSE)
cat("\n\nSample-level data saved to: samples_with_environmental_binary.csv\n")

# Summary by population
cat("\n\n=== Environmental classification summary by population ===\n")
pop_classification <- sample_environmental %>%
  group_by(pop) %>%
  summarise(
    n_samples = n(),
    temp_class = first(mean_annual_temp_binary),
    rainfall_class = first(annual_rainfall_binary),
    wind_class = first(mean_wind_speed_binary),
    .groups = 'drop'
  ) %>%
  mutate(
    temp_label = ifelse(temp_class == 1, "High", "Low"),
    rainfall_label = ifelse(rainfall_class == 1, "High", "Low"),
    wind_label = ifelse(wind_class == 1, "High", "Low")
  )

print(pop_classification)

cat("\n\n=== Cross-tabulation of environmental conditions ===\n")
cat("\nTemperature vs Rainfall:\n")
print(table(
  Temperature = pop_classification$temp_label,
  Rainfall = pop_classification$rainfall_label
))

cat("\nTemperature vs Wind:\n")
print(table(
  Temperature = pop_classification$temp_label,
  Wind = pop_classification$wind_label
))

cat("\n\n=== Analysis complete ===\n")
cat("Files created:\n")
cat("  - transformation_report_median.csv: Median-split parameters\n")
cat("  - transformation_report_threshold.csv: Threshold parameters\n")
cat("  - population_environmental_binary.csv: Population-level binary classifications\n")
cat("  - samples_with_environmental_binary.csv: Sample-level classifications\n")
cat("\nThese binary environmental variables (0/1) can now be used in:\n")
cat("  - Association studies (genotype-environment associations)\n")
cat("  - Population structure analysis\n")
cat("  - Environmental adaptation studies\n")
cat("  - ANOVA or regression models\n")
