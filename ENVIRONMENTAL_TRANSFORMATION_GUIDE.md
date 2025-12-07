# Environmental Variable Transformation Guide

## Overview

This guide explains how to transform continuous environmental variables (temperature, rainfall, wind speed, etc.) into binary values (0 and 1) for use in population genetics and environmental adaptation studies.

## Background

### Why Transform to Binary?

In population genetics studies examining environmental adaptation (such as the paper "Temperature, rainfall and wind variables underlie environmental adaptation in natural populations of Drosophila melanogaster"), continuous environmental variables are often transformed to binary categories for several reasons:

1. **Simplification**: Reduces complex continuous data to interpretable categories
2. **Statistical Power**: Increases power to detect associations in some analyses
3. **Biological Interpretation**: Creates meaningful categories (e.g., "high" vs "low" temperature)
4. **Categorical Analysis**: Enables use of categorical statistical methods (chi-square, logistic regression)
5. **Gene-Environment Associations**: Facilitates testing for genetic variants associated with environmental extremes

## Transformation Methods

### 1. Median-Split Method

**Description**: Divides populations into two equal groups based on the median value.

**Rule**: 
- Values > median = 1 (High)
- Values ≤ median = 0 (Low)

**When to Use**:
- When you want equal-sized groups
- When there's no specific biological threshold
- For exploratory analyses

**Example**:
```r
source("environmental_variable_transformation.R")

result <- transform_environmental_variables(
  data = your_data,
  variables = c("temperature", "rainfall", "wind_speed"),
  method = "median_split"
)
```

**Advantages**:
- Simple and objective
- Always produces balanced groups
- No need for prior knowledge of thresholds

**Disadvantages**:
- Threshold varies by dataset
- May not reflect biological significance
- Treats values near median as categorically different

### 2. Threshold Method

**Description**: Uses predefined, biologically meaningful thresholds.

**Rule**:
- Values > threshold = 1 (High)
- Values ≤ threshold = 0 (Low)

**When to Use**:
- When biological/ecological thresholds are known
- For comparison across studies
- When thresholds have functional significance

**Example**:
```r
result <- transform_environmental_variables(
  data = your_data,
  variables = c("temperature", "rainfall"),
  method = "threshold",
  thresholds = list(
    temperature = 20,   # 20°C threshold
    rainfall = 1000     # 1000mm threshold
  )
)
```

**Advantages**:
- Biologically interpretable
- Consistent across studies
- Can reflect known ecological thresholds

**Disadvantages**:
- Requires domain knowledge
- May produce unbalanced groups
- Threshold choice affects results

### 3. Quartile Method

**Description**: Focuses on environmental extremes by comparing upper and lower quartiles.

**Rule**:
- Values ≥ 75th percentile = 1 (High)
- Values ≤ 25th percentile = 0 (Low)
- Middle 50% = NA or excluded

**When to Use**:
- Focus on extreme environmental conditions
- When intermediate values are less informative
- For detecting adaptation to environmental extremes

**Example**:
```r
result <- transform_environmental_variables(
  data = your_data,
  variables = c("temperature", "rainfall", "wind_speed"),
  method = "quartile"
)
```

**Advantages**:
- Focuses on clear environmental differences
- Reduces noise from intermediate values
- Good for detecting extreme adaptations

**Disadvantages**:
- Excludes middle 50% of data
- Reduces sample size
- May miss linear relationships

### 4. Z-Score Method

**Description**: Based on standardized deviations from the mean.

**Rule**:
- Values above mean (z > 0) = 1 (High)
- Values at or below mean (z ≤ 0) = 0 (Low)

**When to Use**:
- When standardizing across multiple variables
- For combining with other statistical methods
- When mean is biologically relevant

**Example**:
```r
result <- transform_environmental_variables(
  data = your_data,
  variables = c("temperature", "rainfall", "wind_speed"),
  method = "zscore"
)
```

## Practical Workflow

### Step 1: Prepare Your Data

Ensure your data frame contains:
- Population/sample identifiers
- Geographic coordinates (lat, lon)
- Continuous environmental variables

```r
# Example data structure
#   pop    lat       lon      temperature  rainfall  wind_speed
#   24   -33.72   -54.26       18.5        950       14.2
#   28   -31.20   -55.65       19.8        1050      15.8
```

### Step 2: Choose Transformation Method

Consider:
- Research question
- Available biological knowledge
- Sample size
- Distribution of variables

### Step 3: Apply Transformation

```r
source("environmental_variable_transformation.R")

# Transform variables
result <- transform_environmental_variables(
  data = environmental_data,
  variables = c("temperature", "rainfall", "wind_speed"),
  method = "median_split",  # or "threshold", "quartile", "zscore"
  suffix = "_binary"
)

# Get transformed data
transformed_data <- result$data

# Get transformation parameters
transformation_info <- result$transformation_info
```

### Step 4: Document Transformation

```r
# Generate and save transformation report
report <- generate_transformation_report(
  transformation_info,
  output_file = "transformation_report.csv"
)
```

### Step 5: Analyze Binary Variables

Use transformed variables in:
- Association mapping
- ANOVA/regression
- Population differentiation analysis
- Gene-environment correlation studies

## Example Application: Bromus auleticus

The repository includes an example analysis:

```bash
Rscript example_environmental_analysis.R
```

This script:
1. Loads population metadata
2. Simulates environmental data (replace with real climate data)
3. Transforms variables using multiple methods
4. Creates population and sample-level classifications
5. Generates summary statistics and reports

## Extracting Real Environmental Data

For actual analyses, extract environmental variables from climate databases:

### WorldClim (https://www.worldclim.org/)
```r
library(raster)
library(sp)

# Download bioclimatic variables
bio <- getData('worldclim', var='bio', res=2.5)

# Extract for your coordinates
coords <- data.frame(lon = metadata$lon, lat = metadata$lat)
env_values <- extract(bio, coords)
```

### CHELSA (https://chelsa-climate.org/)
- High-resolution climate data
- Temperature, precipitation, and other variables

### Other Sources
- Wind data: Global Wind Atlas
- Historical data: CRU, ERA5
- Local weather stations

## Best Practices

1. **Document Your Method**: Always record which transformation method you used and why

2. **Report Parameters**: Include thresholds, medians, or other parameters in your methods

3. **Check Distributions**: Examine the distribution of binary variables (balance between 0 and 1)

4. **Sensitivity Analysis**: Test multiple methods to ensure results are robust

5. **Biological Relevance**: Consider whether your categorization makes biological sense

6. **Cross-Validation**: If possible, validate thresholds using independent data

7. **Multiple Testing**: Account for multiple testing when using multiple environmental variables

## Citation

When using this methodology, cite relevant papers such as:

- The original Drosophila melanogaster environmental adaptation study
- Papers establishing thresholds for your study system
- Climate data sources used

## Files in This Repository

- `environmental_variable_transformation.R`: Core transformation functions
- `example_environmental_analysis.R`: Example workflow with Bromus auleticus data
- `ENVIRONMENTAL_TRANSFORMATION_GUIDE.md`: This documentation

## Troubleshooting

### Issue: Unbalanced Groups

**Solution**: 
- Use median-split method for balanced groups
- Or adjust threshold to achieve desired balance

### Issue: All Values in One Category

**Solution**: 
- Check if threshold is appropriate for your data range
- Consider using a different method
- Examine data for outliers or errors

### Issue: Missing Values

**Solution**: 
- Functions handle NA values automatically
- Decide whether to exclude or impute missing data
- Document your approach

## Further Reading

1. Population genetics textbooks on quantitative traits
2. Environmental adaptation study designs
3. Statistical methods for gene-environment interactions
4. Papers on environmental niche modeling

## Support

For questions or issues with this implementation, please refer to the repository documentation or contact the maintainers.
