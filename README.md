# SDM – *Alouatta belzebul*

![R](https://img.shields.io/badge/r-%23276DC3.svg?style=flat&logo=r&logoColor=white) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24426806.svg)](https://doi.org/10.5281/zenodo.24426806)

R script for species distribution modelling (SDM) of *Alouatta belzebul* using
the biomod2 package (v4.2-5).

## What the script does
1. Loads occurrence records (literature and field sampling), the study area and
   the environmental rasters
2. Selects predictors using Pearson correlation and VIF
3. Prepares presence and background data
4. Fits and evaluates single models (GLM, RF, ANN, MARS, FDA, CTA, SRE, XGBOOST)
5. Builds a TSS-weighted ensemble and extracts variable importance and response curves
6. Projects the ensemble onto current conditions

## Folder structure
- `1_OCCURRENCES/`: species occurrence records (field and literature)
- `2_ENV_VARIABLES/`: environmental rasters (.tif)
- `3_STUDY_AREA/`: study area shapefile
- `4_RESULTS/`: modelling outputs

## Reproducibility
The script checks and installs the package versions available on the CRAN
snapshot of 2024-06-01 (Posit Public Package Manager).
