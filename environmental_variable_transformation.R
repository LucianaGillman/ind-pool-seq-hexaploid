#!/usr/bin/env Rscript
# Environmental Variable Transformation Script
# 
# This script implements the methodology for transforming continuous environmental
# variables (temperature, rainfall, wind) into binary values (0 and 1).
#
# Methodology based on common practices in population genetics studies, including
# the approach from "Temperature, rainfall and wind variables underlie environmental
# adaptation in natural populations of Drosophila melanogaster"
#
# Binary transformation methods:
# 1. Median-split: Values above median = 1, below median = 0
# 2. Threshold-based: Values above defined threshold = 1, below = 0
# 3. Quartile-based: Upper quartile = 1, lower quartile = 0

# Load required libraries
library(dplyr)
library(tidyr)

# Function to perform median-split binary transformation
# This is the most common method used in environmental adaptation studies
median_split_transform <- function(variable, na_value = NA) {
  # Calculate median, removing NA values
  med <- median(variable, na.rm = TRUE)
  
  # Transform to binary: 1 if above median, 0 if below or equal
  binary <- ifelse(is.na(variable), na_value,
                   ifelse(variable > med, 1, 0))
  
  return(list(
    binary = binary,
    threshold = med,
    method = "median_split"
  ))
}

# Function to perform threshold-based transformation
threshold_transform <- function(variable, threshold, na_value = NA) {
  # Transform to binary based on provided threshold
  binary <- ifelse(is.na(variable), na_value,
                   ifelse(variable > threshold, 1, 0))
  
  return(list(
    binary = binary,
    threshold = threshold,
    method = "threshold"
  ))
}

# Function to perform quartile-based transformation
# Useful for focusing on extreme environmental conditions
quartile_transform <- function(variable, na_value = NA) {
  # Calculate quartiles
  q1 <- quantile(variable, 0.25, na.rm = TRUE)
  q3 <- quantile(variable, 0.75, na.rm = TRUE)
  
  # Transform: upper quartile = 1, lower quartile = 0, middle = NA or 0.5
  binary <- ifelse(is.na(variable), na_value,
                   ifelse(variable >= q3, 1,
                          ifelse(variable <= q1, 0, NA)))
  
  return(list(
    binary = binary,
    q1 = q1,
    q3 = q3,
    method = "quartile"
  ))
}

# Function to perform z-score based transformation
# Values above mean (z > 0) = 1, below mean (z <= 0) = 0
zscore_transform <- function(variable, na_value = NA) {
  # Calculate z-scores
  # Note: scale() returns a matrix with one column; we convert to vector
  # for compatibility with ifelse()
  z <- scale(variable, center = TRUE, scale = TRUE)
  z <- as.vector(z)  # Convert matrix to vector
  
  # Transform based on z-score
  binary <- ifelse(is.na(variable), na_value,
                   ifelse(z > 0, 1, 0))
  
  mean_val <- mean(variable, na.rm = TRUE)
  sd_val <- sd(variable, na.rm = TRUE)
  
  return(list(
    binary = binary,
    mean = mean_val,
    sd = sd_val,
    method = "zscore"
  ))
}

# Main function to transform environmental variables
transform_environmental_variables <- function(data, 
                                             variables,
                                             method = "median_split",
                                             thresholds = NULL,
                                             suffix = "_binary") {
  #' Transform environmental variables to binary (0/1)
  #' 
  #' @param data Data frame containing environmental variables
  #' @param variables Character vector of variable names to transform
  #' @param method Transformation method: "median_split", "threshold", 
  #'               "quartile", or "zscore"
  #' @param thresholds Named list of thresholds (for threshold method)
  #' @param suffix Suffix to add to binary variable names
  #' 
  #' @return Data frame with original data and binary transformed variables
  
  result_data <- data
  transformation_info <- list()
  
  for (var in variables) {
    if (!var %in% names(data)) {
      warning(paste("Variable", var, "not found in data"))
      next
    }
    
    # Apply transformation based on method
    if (method == "median_split") {
      trans <- median_split_transform(data[[var]])
    } else if (method == "threshold") {
      if (is.null(thresholds) || !var %in% names(thresholds)) {
        warning(paste("No threshold specified for", var, "- using median"))
        trans <- median_split_transform(data[[var]])
      } else {
        trans <- threshold_transform(data[[var]], thresholds[[var]])
      }
    } else if (method == "quartile") {
      trans <- quartile_transform(data[[var]])
    } else if (method == "zscore") {
      trans <- zscore_transform(data[[var]])
    } else {
      stop(paste("Unknown method:", method))
    }
    
    # Add binary variable to data
    new_var_name <- paste0(var, suffix)
    result_data[[new_var_name]] <- trans$binary
    
    # Store transformation information
    transformation_info[[var]] <- trans[names(trans) != "binary"]
  }
  
  return(list(
    data = result_data,
    transformation_info = transformation_info
  ))
}

# Function to generate transformation report
generate_transformation_report <- function(transformation_info, output_file = NULL) {
  #' Generate a report of transformation parameters
  #' 
  #' @param transformation_info List of transformation information
  #' @param output_file Optional file path to save report
  #' 
  #' @return Data frame with transformation parameters
  
  report_rows <- list()
  
  for (var_name in names(transformation_info)) {
    info <- transformation_info[[var_name]]
    
    row <- data.frame(
      Variable = var_name,
      Method = info$method,
      stringsAsFactors = FALSE
    )
    
    # Add method-specific information
    if (info$method == "median_split" || info$method == "threshold") {
      row$Threshold <- info$threshold
    } else if (info$method == "quartile") {
      row$Q1 <- info$q1
      row$Q3 <- info$q3
    } else if (info$method == "zscore") {
      row$Mean <- info$mean
      row$SD <- info$sd
    }
    
    report_rows[[var_name]] <- row
  }
  
  report <- bind_rows(report_rows)
  
  if (!is.null(output_file)) {
    write.csv(report, output_file, row.names = FALSE)
    cat("Transformation report saved to:", output_file, "\n")
  }
  
  return(report)
}

# Example usage function
example_usage <- function() {
  #' Example demonstrating how to use the transformation functions
  
  # Create sample environmental data
  set.seed(123)
  n <- 100
  sample_data <- data.frame(
    population = paste0("Pop", 1:n),
    temperature = rnorm(n, mean = 20, sd = 5),
    rainfall = rnorm(n, mean = 800, sd = 200),
    wind_speed = rnorm(n, mean = 15, sd = 4)
  )
  
  cat("=== Example 1: Median-split transformation ===\n")
  result1 <- transform_environmental_variables(
    data = sample_data,
    variables = c("temperature", "rainfall", "wind_speed"),
    method = "median_split"
  )
  print(head(result1$data))
  print(generate_transformation_report(result1$transformation_info))
  
  cat("\n=== Example 2: Threshold-based transformation ===\n")
  result2 <- transform_environmental_variables(
    data = sample_data,
    variables = c("temperature", "rainfall", "wind_speed"),
    method = "threshold",
    thresholds = list(temperature = 20, rainfall = 750, wind_speed = 15)
  )
  print(head(result2$data))
  print(generate_transformation_report(result2$transformation_info))
  
  cat("\n=== Example 3: Quartile-based transformation ===\n")
  result3 <- transform_environmental_variables(
    data = sample_data,
    variables = c("temperature", "rainfall", "wind_speed"),
    method = "quartile"
  )
  print(head(result3$data))
  print(generate_transformation_report(result3$transformation_info))
}

# Main execution when script is run directly
if (!interactive()) {
  cat("Environmental Variable Transformation Script\n")
  cat("=============================================\n\n")
  cat("This script provides functions to transform continuous environmental\n")
  cat("variables into binary (0/1) values using various methods.\n\n")
  cat("Run example_usage() to see demonstration.\n")
}
