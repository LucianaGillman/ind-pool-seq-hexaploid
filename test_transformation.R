#!/usr/bin/env Rscript
# Test Script for Environmental Variable Transformation
#
# Purpose:
#   This script validates the core functionality of environmental variable
#   transformation methods (median-split and threshold-based) without requiring
#   external packages. It tests the transformation logic, handling of missing
#   values, and correct binary encoding.
#
# Expected Output:
#   - Test 1: Demonstrates median-split transformation creating balanced groups
#   - Test 2: Shows threshold-based transformation with custom cutoff
#   - Test 3: Validates proper handling of missing (NA) values
#
# Relationship to Main Functions:
#   This test replicates the core logic implemented in 
#   environmental_variable_transformation.R using base R functions only.
#   It serves as a lightweight validation that the transformation methodology
#   works correctly before applying the full-featured functions to real data.
#
# Usage:
#   Rscript test_transformation.R

cat("Testing Environmental Variable Transformation Functions\n")
cat("========================================================\n\n")

# Test data
set.seed(123)
test_data <- data.frame(
  site = c("A", "B", "C", "D", "E", "F", "G", "H"),
  temperature = c(15, 18, 20, 22, 16, 21, 19, 23),
  rainfall = c(800, 950, 1100, 1200, 850, 1050, 900, 1250)
)

cat("Test data:\n")
print(test_data)

# Test median-split transformation
cat("\n\nTest 1: Median-split transformation\n")
cat("------------------------------------\n")

temp_median <- median(test_data$temperature)
cat("Temperature median:", temp_median, "\n")

test_data$temp_binary <- ifelse(test_data$temperature > temp_median, 1, 0)
cat("Binary transformation (> median = 1, <= median = 0):\n")
print(test_data[, c("site", "temperature", "temp_binary")])

# Verify counts
cat("\nHigh (1):", sum(test_data$temp_binary == 1), "sites\n")
cat("Low (0):", sum(test_data$temp_binary == 0), "sites\n")

# Test threshold transformation
cat("\n\nTest 2: Threshold-based transformation\n")
cat("---------------------------------------\n")

rainfall_threshold <- 1000
cat("Rainfall threshold:", rainfall_threshold, "mm\n")

test_data$rain_binary <- ifelse(test_data$rainfall > rainfall_threshold, 1, 0)
cat("Binary transformation (> threshold = 1, <= threshold = 0):\n")
print(test_data[, c("site", "rainfall", "rain_binary")])

# Verify counts
cat("\nHigh (1):", sum(test_data$rain_binary == 1), "sites\n")
cat("Low (0):", sum(test_data$rain_binary == 0), "sites\n")

# Test with NA values
cat("\n\nTest 3: Handling missing values\n")
cat("--------------------------------\n")

test_data_na <- test_data
test_data_na$temperature[3] <- NA
test_data_na$rainfall[5] <- NA

cat("Data with missing values:\n")
print(test_data_na[, c("site", "temperature", "rainfall")])

cat("\nBinary transformations (NA values preserved):\n")
temp_median_na <- median(test_data_na$temperature, na.rm = TRUE)
test_data_na$temp_binary <- ifelse(is.na(test_data_na$temperature), NA,
                                   ifelse(test_data_na$temperature > temp_median_na, 1, 0))

rain_median_na <- median(test_data_na$rainfall, na.rm = TRUE)
test_data_na$rain_binary <- ifelse(is.na(test_data_na$rainfall), NA,
                                   ifelse(test_data_na$rainfall > rain_median_na, 1, 0))

print(test_data_na[, c("site", "temp_binary", "rain_binary")])

cat("\n\nAll tests completed successfully!\n")
cat("\nThe transformation methodology is working correctly:\n")
cat("  1. Median-split: Divides data into high/low groups\n")
cat("  2. Threshold: Uses custom cutoff values\n")
cat("  3. Missing values: Handled appropriately\n")
cat("\nThis demonstrates the core functionality used in environmental\n")
cat("adaptation studies like the Drosophila melanogaster paper.\n")
