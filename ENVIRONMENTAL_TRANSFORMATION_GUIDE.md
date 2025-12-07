# Environmental Variable Transformation Guide

## Overview

This guide explains how environmental variables (eco-regions) are transformed into binary (0/1) categories for statistical analysis in the Bromus auleticus genomic diversity study.

## Background

The repository contains:
- **Spatial data**: Eco-region shapefile (`eco-regiones_27-4-2012/`) with 8 eco-regions of Uruguay
- **Sample coordinates**: Latitude and longitude for 296 samples in `metadata.csv`
- **5 populations**: Samples grouped into populations (24, 28, 50, 87, 88)

## The Problem

Environmental variables like eco-regions are categorical data. To use them in statistical analyses (ANOVA, regression, etc.), they need to be transformed into binary variables.

## The Solution: Spatial Overlay Analysis

### Step 1: Spatial Join
Use the sample coordinates to determine which eco-region each sample belongs to by performing a spatial join (point-in-polygon operation).

### Step 2: Binary Transformation
Transform the categorical eco-region assignment into binary (0/1) variables using one of these methods:

#### Method 1: One-Hot Encoding
Create one binary column per eco-region:
- Value = 1 if sample is in that eco-region
- Value = 0 if sample is NOT in that eco-region

Example:
```
Sample | ECO_REGION               | eco_Graven_Laguna | eco_Escudo | eco_Sedimentaria
-------|--------------------------|-------------------|------------|------------------
001    | Graven de la Laguna Merín|        1          |     0      |        0
002    | Escudo Cristalino        |        0          |     1      |        0
```

#### Method 2: Grouped Categories
Combine similar eco-regions into broader binary categories:

1. **Coastal vs Inland**
   - Coastal (1): Graven de la Laguna Merín, Sierras del Este
   - Inland (0): All others

2. **Sedimentary Basin vs Other**
   - Sedimentary (1): Cuenca sedimentaria Gondwánica, Cuenca sedimentaria del Oeste
   - Other (0): All others

3. **Highland vs Lowland**
   - Highland (1): Cuesta Basáltica, Sierras del Este, Escudo Cristalino
   - Lowland (0): All others

### Step 3: Aggregation by Population
Since multiple samples belong to each population, calculate:
- Proportion of samples in each eco-region per population
- Most common eco-region for each population

## Eco-regions in the Dataset

The shapefile contains 8 eco-regions of Uruguay:

1. **Cuesta Basáltica** - Basaltic ridge formation
2. **Cuenca sedimentaria Gondwánica** - Gondwanic sedimentary basin
3. **Sierras del Este** - Eastern hills
4. **Cuenca sedimentaria del Oeste** - Western sedimentary basin
5. **Graven de la Laguna Merín** - Merín Lagoon graben
6. **Graven del Santa Lucía** - Santa Lucía graben
7. **Escudo Cristalino** - Crystalline shield

## Distribution in This Study

From the analysis, samples are distributed across 4 eco-regions:

| Population | Eco-region                      | Number of Samples |
|------------|---------------------------------|-------------------|
| 24         | Graven de la Laguna Merín       | 60                |
| 28         | Cuenca sedimentaria Gondwánica  | 60                |
| 50         | Escudo Cristalino               | 58                |
| 87         | Cuenca sedimentaria Gondwánica  | 58                |
| 88         | Graven del Santa Lucía          | 60                |

**Key finding**: Each population is entirely within a single eco-region, making the binary variables perfectly correlated with population identity in this dataset.

## Applications

Binary environmental variables can be used for:

1. **ANOVA**: Test if genetic diversity differs between eco-regions
2. **Regression**: Model genetic diversity as a function of environmental categories
3. **PCA/MDS**: Include environmental variables as supplementary data
4. **RDA (Redundancy Analysis)**: Relate genetic variation to environmental predictors
5. **Mantel tests**: Test correlation between genetic and environmental distance matrices

## Continuous Variables Alternative

Instead of eco-regions, you can extract continuous environmental variables (temperature, precipitation, elevation) at each sample location and then:

1. **Median split**: High (1) vs Low (0) based on median value
2. **Quartile split**: Top quartile (1) vs Bottom quartile (0)
3. **Threshold-based**: Above/below ecologically meaningful threshold

## Implementation

Two scripts are provided to perform this transformation:

- `environmental_variable_transformation.py` (Python)
- `environmental_variable_transformation.R` (R)

Both scripts:
- Perform spatial overlay
- Create binary variables
- Generate summaries
- Export results to CSV
- Create visualizations

See the main README.md for usage instructions.

## References

This methodology follows standard practices in landscape genetics and ecological genetics:

- Spatial overlay: Common GIS operation for assigning environmental data to sample locations
- Binary transformation: Standard approach for incorporating categorical predictors in statistical models
- One-hot encoding: Widely used in machine learning and statistics for categorical variables

## Notes

- The original shapefile uses "Graven" (Spanish spelling) instead of "Graben" (German/English geological term)
- All accented characters in eco-region names are transliterated in column names (á→a, é→e, etc.)
- Missing data is handled gracefully (samples outside all polygons get NA values)
