This repository contains all scripts used for the paper "Optimizing Genomic Diversity Assessments for Conservation of Bromus auleticus Using Individual and Pooled Sequencing". These scripts allow comparison of individual sequencing (ind-seq) and pooled sequencing (pool-seq) techniques.

## Environmental Variable Transformation

### How to Transform Environmental Variables to Binary (0/1) Categories

The repository includes eco-region spatial data (`eco-regiones_27-4-2012/`) and sample coordinates (`metadata.csv`). To transform environmental variables (eco-regions) into binary categories, two equivalent scripts are provided (Python and R versions).

These scripts demonstrate the complete workflow:

1. **Spatial Overlay**: Uses the sample coordinates (latitude/longitude) to determine which eco-region each sample belongs to by performing a spatial join with the eco-region shapefile.

2. **Binary Transformation Methods**:
   - **One-hot encoding**: Creates one binary column per eco-region (1 = sample in that region, 0 = not in that region)
   - **Grouped categories**: Combines similar eco-regions into broader binary categories based on characteristics:
     - Coastal vs Inland
     - Sedimentary basin vs Other geological formations
     - Highland vs Lowland regions

3. **Population-level Summary**: Aggregates binary variables by population since multiple samples belong to each population.

### Available Scripts

Two equivalent scripts are provided:

#### Python Version
**File:** `environmental_variable_transformation.py`

Required packages:
```bash
pip install geopandas matplotlib pandas
```

Usage:
```bash
python environmental_variable_transformation.py
```

#### R Version
**File:** `environmental_variable_transformation.R`

Required packages:
```r
install.packages(c("sf", "dplyr", "ggplot2"))
```

Usage:
```r
source("environmental_variable_transformation.R")
```

### Output Files

- `metadata_with_binary_ecoregions.csv`: Sample-level data with binary eco-region variables
- `population_ecoregion_summary.csv`: Population-level summary of eco-region membership
- `sample_distribution_by_ecoregion.png`: Visualization of sample distribution

### Eco-regions in the Dataset

The shapefile contains 8 eco-regions of Uruguay:
1. Cuesta Basáltica
2. Cuenca sedimentaria Gondwánica
3. Sierras del Este
4. Cuenca sedimentaria del Oeste
5. Graven de la Laguna Merín
6. Graven del Santa Lucía
7. Escudo Cristalino

Each can be transformed into a binary variable for statistical analysis.
