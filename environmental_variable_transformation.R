# Environmental Variable Transformation to Binary Categories
# This script demonstrates how to transform eco-region environmental data 
# to binary (0/1) categories based on spatial overlay analysis

# Load required libraries
library(sf)           # For spatial data handling
library(dplyr)        # For data manipulation

# Set working directory (adjust as needed)
# setwd("/path/to/your/working/directory")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Load metadata with sample coordinates
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# Load eco-region shapefile
eco_regions <- st_read("eco-regiones_27-4-2012/eco-regiones_27-4-2012.shp")

# ============================================================================
# 2. PREPARE SPATIAL DATA
# ============================================================================

# Convert metadata to spatial points object
# Assuming coordinates are in WGS84 (EPSG:4326)
metadata_sf <- st_as_sf(metadata, 
                        coords = c("lon", "lat"), 
                        crs = 4326)

# Ensure both datasets have the same CRS
metadata_sf <- st_transform(metadata_sf, st_crs(eco_regions))

# ============================================================================
# 3. SPATIAL OVERLAY: Assign eco-region to each sample
# ============================================================================

# Perform spatial join to determine which eco-region each sample belongs to
metadata_with_ecoregion <- st_join(metadata_sf, eco_regions)

# Convert back to regular dataframe
metadata_df <- as.data.frame(metadata_with_ecoregion)
metadata_df <- metadata_df %>% select(-geometry)

# ============================================================================
# 4. TRANSFORM ECO-REGIONS TO BINARY (0/1) CATEGORIES
# ============================================================================

# Method 1: Create binary variable for each eco-region
# Each column represents whether a sample belongs to that eco-region (1) or not (0)

# Get unique eco-region names
unique_ecoregions <- unique(eco_regions$ECO_REGION)

# Create binary columns for each eco-region
for (region in unique_ecoregions) {
  # Create column name (clean special characters)
  col_name <- gsub(" ", "_", region)
  # Use chartr for efficient character replacement
  col_name <- chartr("áéíóú", "aeiou", col_name)
  col_name <- paste0("eco_", col_name)
  
  # Create binary variable: 1 if sample is in this eco-region, 0 otherwise
  metadata_df[[col_name]] <- ifelse(metadata_df$ECO_REGION == region, 1, 0)
}

# ============================================================================
# 5. METHOD 2: Group eco-regions into broader categories
# ============================================================================

# Example: Create binary variables based on eco-region characteristics
# You can customize this based on your research questions

# Example 1: Coastal vs Inland
metadata_df$coastal <- ifelse(metadata_df$ECO_REGION %in% 
                               c("Graven de la Laguna Merín", "Sierras del Este"), 
                               1, 0)

# Example 2: Sedimentary basin vs Other
metadata_df$sedimentary_basin <- ifelse(grepl("sedimentaria", 
                                              metadata_df$ECO_REGION, 
                                              ignore.case = TRUE), 
                                        1, 0)

# Example 3: High vs Low elevation (based on eco-region names)
metadata_df$highland <- ifelse(metadata_df$ECO_REGION %in% 
                               c("Cuesta Basáltica", "Sierras del Este", "Escudo Cristalino"), 
                               1, 0)

# ============================================================================
# 6. SUMMARIZE BY POPULATION
# ============================================================================

# Since multiple samples belong to the same population (pop column),
# we can summarize the binary variables by population

# Calculate the proportion of samples in each eco-region per population
pop_summary <- metadata_df %>%
  group_by(pop) %>%
  summarise(
    n_samples = n(),
    eco_region = first(ECO_REGION),  # Most common eco-region
    across(starts_with("eco_"), mean),  # Proportion in each eco-region
    coastal = mean(coastal),
    sedimentary_basin = mean(sedimentary_basin),
    highland = mean(highland)
  )

# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

# Save the transformed data with binary variables
write.csv(metadata_df, "metadata_with_binary_ecoregions.csv", row.names = FALSE)

# Save population-level summary
write.csv(pop_summary, "population_ecoregion_summary.csv", row.names = FALSE)

# Print summary
cat("\n=== Eco-region Distribution ===\n")
print(table(metadata_df$ECO_REGION))

cat("\n=== Population Summary ===\n")
print(pop_summary)

# ============================================================================
# 8. OPTIONAL: Create visualization
# ============================================================================

# Plot samples colored by eco-region
if (require(ggplot2)) {
  p <- ggplot(metadata_df, aes(x = lon, y = lat, color = ECO_REGION)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "Sample Distribution by Eco-region",
         x = "Longitude", y = "Latitude")
  
  ggsave("sample_distribution_by_ecoregion.png", p, width = 10, height = 8)
  print(p)
}

# ============================================================================
# NOTES:
# ============================================================================
# 
# The binary transformation can be done in multiple ways depending on your
# research question:
# 
# 1. One-hot encoding: Create one binary column per eco-region (done above)
# 
# 2. Group similar eco-regions: Combine eco-regions based on:
#    - Ecological characteristics (e.g., coastal vs inland)
#    - Geological features (e.g., sedimentary vs crystalline)
#    - Climate zones
#    - Elevation
# 
# 3. Environmental variables: Instead of eco-regions, you can extract 
#    continuous environmental variables (temperature, precipitation, etc.)
#    and then discretize them into binary categories using:
#    - Median split (above/below median = 1/0)
#    - Quartile split (high vs low quartiles = 1/0)
#    - Ecologically meaningful thresholds
# 
# 4. For statistical analysis (e.g., ANOVA, regression), you might want to:
#    - Use the binary variables as predictors
#    - Test association between genetic diversity and eco-region
#    - Compare populations from different eco-regions
# 
# ============================================================================
