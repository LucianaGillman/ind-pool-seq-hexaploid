This repository contains all scripts used for the paper "Optimizing Genomic Diversity Assessments for Conservation of Bromus auleticus Using Individual and Pooled Sequencing". These scripts allow comparison of individual sequencing (ind-seq) and pooled sequencing (pool-seq) techniques.

## Environmental Variable Transformation

This repository now includes tools for transforming continuous environmental variables (temperature, rainfall, wind speed) into binary values (0 and 1), following methodologies used in population genetics and environmental adaptation studies.

### Quick Start

```r
# Source the transformation functions
source("environmental_variable_transformation.R")

# Run the example analysis
source("example_environmental_analysis.R")
```

### Key Features

- **Multiple Transformation Methods**:
  - Median-split: Divides populations into equal high/low groups
  - Threshold-based: Uses biologically meaningful cutoffs
  - Quartile-based: Focuses on environmental extremes
  - Z-score based: Standardized transformation

- **Comprehensive Documentation**: See `ENVIRONMENTAL_TRANSFORMATION_GUIDE.md` for detailed explanations

- **Example Workflow**: `example_environmental_analysis.R` demonstrates the complete process

### Files

- `environmental_variable_transformation.R`: Core transformation functions
- `example_environmental_analysis.R`: Example using Bromus auleticus data
- `ENVIRONMENTAL_TRANSFORMATION_GUIDE.md`: Complete methodology guide

### Use Cases

- Gene-environment association studies
- Environmental adaptation analysis
- Population structure analysis with environmental context
- Genotype-phenotype-environment relationships
