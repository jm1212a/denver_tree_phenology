# Urban Tree Phenology Analysis - Denver, Colorado

Analysis of climate impacts on urban tree phenology across 12 deciduous species in Denver, using satellite imagery and XGBoost modeling to predict four phenological stages (2018-2024).

**Author:** Jared Marint 
**Institution:** American University  
**Program:** M.S. Data Science  
**Completion:** December 2025

## Project Overview

This capstone project analyzes how climate variables influence the timing of seasonal changes (phenology) in urban trees. Using PlanetScope satellite imagery, NOAA climate data, and Denver's tree inventory, the study models phenological transitions for 3,450 trees across 12 deciduous species over a 7-year period.

### Interactive Shiny Application
An interactive Shiny application (`denver_trees_app/`) was developed to facilitate sample selection and quality control during the research process.

## Repository Structure
```
Denver_Tree_Phenology/
├── Final_Draft_Species_S...n_Denver_Colorado.qmd    # Main analysis document
├── Final_Draft--_Species_S...n_Denver-_Colorado.pdf # Rendered output
├── denver_tree_phenology.Rproj                      # R Project file
├── README.md                                        
│
├── analysis/
│   ├── modeling/          # XGBoost model training and evaluation
│   ├── graphs/            # Visualization scripts
│
├── data/
│   └── processed_trees_ndvi_lc_climate.rds          # Final processed dataset (1.5 GB)
│
├── data_raw/
│   ├── noaa_daily_weather_obs/    # NOAA climate data
│   ├── landcover/                  # Land cover data
│   ├── ndvi/                       # NDVI time series
│   ├── denver_tree/                # Denver tree inventory
│   ├── tree_leaf_taxanomy/         # Species characteristics
│   └── [other raw data folders]
│
├── denver_trees_app/      # Interactive Shiny application for QC and exploration
│
├── tables/                # Summary tables for paper
├── images/                # Figures and plots
│
├── literature_review/
│   ├── biblio.bib         # Bibliography
│   └── chicago-author-date.csl    # Citation style
│
└── R/                     # Custom R functions and utilities
```

## Data Sources

1. **Satellite Imagery:** PlanetScope NDVI (2018-2024)
2. **Climate Data:** NOAA daily weather observations (2018-2024) and climate baseline
3. **Tree Inventory:** City and County of Denver Open Data Catalog
4. **Sample:** 3,450 trees, 12 deciduous species, 24,150 tree-year observations

## Methods

**Phenological Stages Analyzed:**
1. Greenup (spring leaf emergence)
2. Maturity (peak greenness)
3. Senescence (fall color change)
4. Dormancy (leaf drop)

**Modeling Approach:**
- XGBoost gradient boosting models
- Species-specific and pooled models
- Climate predictors: temperature, precipitation, growing degree days
- Phenological curve fitting using `phenofit` package
- Parallel processing for large-scale NDVI time series analysis

## Requirements

### R Version
- R ≥ 4.0.0

### Key R Packages
```r
# Data manipulation
library(tidyverse)

# Geospatial analysis
library(sf)
library(terra)
library(leaflet)

# Phenology curvefitting and metric extraction
library(phenofit)

# Machine learning
library(xgboost)
library(caret)

# Data visualtion
library(shiny)

```

## How to Run

1. **Clone the repository**
```bash
   git clone https://github.com/jm1212a/denver_tree_phenology
   cd Denver_Tree_Phenology
```

2. **Open R Project**
```r
   # Open denver_tree_phenology.Rproj in RStudio
```

3. **Install dependencies**
```r
   install.packages(c("tidyverse", "sf", "terra", "xgboost", 
                      "leaflet", "kableExtra", "gt", "phenofit"))
```

4. **Run analysis**
   - Main analysis document: `Final_Draft_Species_S...n_Denver_Colorado.qmd`
   - Render with Quarto to generate PDF output

**Note:** The processed dataset (`processed_trees_ndvi_lc_climate.rds`) is 1.5 GB. Raw data processing scripts are available in `data_raw/` subdirectories.

## Output Files

- **Paper:** `Final_Draft--_Species_S...n_Denver-_Colorado.pdf`
- **Tables:** Located in `tables/` directory
- **Figures:** Located in `images/` directory

## Contact

Jared Martin  
jaredm1880@gamil.com 

## License

MIT License - Feel free to use and adapt this code for research and educational purposes.
