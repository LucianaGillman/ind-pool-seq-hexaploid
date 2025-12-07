#!/usr/bin/env python3
"""
Environmental Variable Transformation to Binary Categories
This script demonstrates how to transform eco-region environmental data 
to binary (0/1) categories based on spatial overlay analysis
"""

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Load metadata with sample coordinates
metadata = pd.read_csv("metadata.csv")

# Load eco-region shapefile
eco_regions = gpd.read_file("eco-regiones_27-4-2012/eco-regiones_27-4-2012.shp")

print("=== Data Loaded ===")
print(f"Number of samples: {len(metadata)}")
print(f"Number of eco-regions: {len(eco_regions)}")
print(f"\nEco-regions in shapefile:")
print(eco_regions['ECO_REGION'].unique())

# ============================================================================
# 2. PREPARE SPATIAL DATA
# ============================================================================

# Convert metadata to GeoDataFrame with Point geometries
geometry = [Point(xy) for xy in zip(metadata['lon'], metadata['lat'])]
metadata_gdf = gpd.GeoDataFrame(metadata, geometry=geometry, crs="EPSG:4326")

# Ensure both datasets have the same CRS
metadata_gdf = metadata_gdf.to_crs(eco_regions.crs)

# ============================================================================
# 3. SPATIAL OVERLAY: Assign eco-region to each sample
# ============================================================================

# Perform spatial join to determine which eco-region each sample belongs to
metadata_with_ecoregion = gpd.sjoin(metadata_gdf, eco_regions, how="left", predicate="within")

# Convert back to regular DataFrame and clean up
metadata_df = pd.DataFrame(metadata_with_ecoregion)
metadata_df = metadata_df.drop(columns=['geometry', 'index_right'], errors='ignore')

print(f"\n=== Eco-region Assignment ===")
print(f"Samples with assigned eco-region: {metadata_df['ECO_REGION'].notna().sum()}")
print(f"\nDistribution of samples by eco-region:")
print(metadata_df['ECO_REGION'].value_counts())

# ============================================================================
# 4. TRANSFORM ECO-REGIONS TO BINARY (0/1) CATEGORIES
# ============================================================================

# Method 1: Create binary variable for each eco-region
# Each column represents whether a sample belongs to that eco-region (1) or not (0)

# Get unique eco-region names (excluding NaN)
unique_ecoregions = metadata_df['ECO_REGION'].dropna().unique()

print("\n=== Creating Binary Variables ===")

# Create binary columns for each eco-region
for region in unique_ecoregions:
    # Create column name (clean special characters)
    col_name = region.replace(" ", "_").replace("á", "a").replace("é", "e")
    col_name = f"eco_{col_name}"
    
    # Create binary variable: 1 if sample is in this eco-region, 0 otherwise
    metadata_df[col_name] = (metadata_df['ECO_REGION'] == region).astype(int)
    print(f"  Created: {col_name} (n={metadata_df[col_name].sum()})")

# ============================================================================
# 5. METHOD 2: Group eco-regions into broader categories
# ============================================================================

print("\n=== Creating Grouped Binary Variables ===")

# Example 1: Coastal vs Inland
coastal_regions = ["Graven de la Laguna Merín", "Sierras del Este"]
metadata_df['coastal'] = metadata_df['ECO_REGION'].isin(coastal_regions).astype(int)
print(f"  Coastal: {metadata_df['coastal'].sum()} samples")

# Example 2: Sedimentary basin vs Other
metadata_df['sedimentary_basin'] = metadata_df['ECO_REGION'].str.contains(
    'sedimentaria', case=False, na=False).astype(int)
print(f"  Sedimentary basin: {metadata_df['sedimentary_basin'].sum()} samples")

# Example 3: High vs Low elevation (based on eco-region names)
highland_regions = ["Cuesta Basáltica", "Sierras del Este", "Escudo Cristalino"]
metadata_df['highland'] = metadata_df['ECO_REGION'].isin(highland_regions).astype(int)
print(f"  Highland: {metadata_df['highland'].sum()} samples")

# ============================================================================
# 6. SUMMARIZE BY POPULATION
# ============================================================================

# Since multiple samples belong to the same population (pop column),
# we can summarize the binary variables by population

# Get all binary eco-region columns
eco_columns = [col for col in metadata_df.columns if col.startswith('eco_')]

# Calculate the proportion of samples in each eco-region per population
pop_summary = metadata_df.groupby('pop').agg({
    'id': 'count',  # Number of samples
    'ECO_REGION': lambda x: x.mode()[0] if len(x.mode()) > 0 else None,  # Most common
    **{col: 'mean' for col in eco_columns},  # Proportion in each eco-region
    'coastal': 'mean',
    'sedimentary_basin': 'mean',
    'highland': 'mean'
}).rename(columns={'id': 'n_samples'})

print("\n=== Population Summary ===")
print(pop_summary)

# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

# Save the transformed data with binary variables
metadata_df.to_csv("metadata_with_binary_ecoregions.csv", index=False)
print("\n✓ Saved: metadata_with_binary_ecoregions.csv")

# Save population-level summary
pop_summary.to_csv("population_ecoregion_summary.csv")
print("✓ Saved: population_ecoregion_summary.csv")

# ============================================================================
# 8. CREATE VISUALIZATION
# ============================================================================

try:
    # Create a simple plot showing sample distribution by eco-region
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot eco-regions
    eco_regions.plot(ax=ax, color='lightgray', edgecolor='black', alpha=0.5)
    
    # Plot samples colored by population
    for pop in metadata_df['pop'].unique():
        pop_data = metadata_df[metadata_df['pop'] == pop]
        ax.scatter(pop_data['lon'], pop_data['lat'], 
                  label=f'Pop {pop}', s=50, alpha=0.7)
    
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Sample Distribution by Population')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sample_distribution_by_population.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: sample_distribution_by_population.png")
    
except Exception as e:
    print(f"Warning: Could not create visualization: {e}")

# ============================================================================
# 9. SUMMARY STATISTICS
# ============================================================================

print("\n" + "="*70)
print("SUMMARY OF BINARY TRANSFORMATION")
print("="*70)

print("\nBinary variables created:")
binary_cols = [col for col in metadata_df.columns if col.startswith('eco_') or 
               col in ['coastal', 'sedimentary_basin', 'highland']]
for col in binary_cols:
    n_ones = metadata_df[col].sum()
    n_zeros = (metadata_df[col] == 0).sum()
    print(f"  {col:40s}: 1s={n_ones:3d}, 0s={n_zeros:3d}")

print("\n" + "="*70)
print("NOTES:")
print("="*70)
print("""
The binary transformation can be done in multiple ways depending on your
research question:

1. One-hot encoding: Create one binary column per eco-region (done above)

2. Group similar eco-regions: Combine eco-regions based on:
   - Ecological characteristics (e.g., coastal vs inland)
   - Geological features (e.g., sedimentary vs crystalline)
   - Climate zones
   - Elevation

3. Environmental variables: Instead of eco-regions, you can extract 
   continuous environmental variables (temperature, precipitation, etc.)
   and then discretize them into binary categories using:
   - Median split (above/below median = 1/0)
   - Quartile split (high vs low quartiles = 1/0)
   - Ecologically meaningful thresholds

4. For statistical analysis (e.g., ANOVA, regression), you might want to:
   - Use the binary variables as predictors
   - Test association between genetic diversity and eco-region
   - Compare populations from different eco-regions
""")

print("\nTransformation complete!")
