# ==============================================================================
# Species distribution modelling of Alouatta belzebul with biomod2
# ==============================================================================
# Workflow:
#   1. Load occurrence records, study area and environmental rasters
#   2. Select predictors (Pearson correlation, then VIF)
#   3. Prepare presence / background data
#   4. Fit single models and evaluate them
#   5. Build the ensemble model, extract importance and response curves
#   6. Project the ensemble onto current conditions
# ==============================================================================


# ==============================================================================
# FOLDER STRUCTURE
# ==============================================================================
# The script runs from a single project folder (the working directory, set in the
# WORKING DIRECTORY section below). This folder must contain four subfolders with
# exactly these names:
#
#   project_folder/
#   |-- 1_OCCURRENCES/    Species occurrence records, both from field sampling and
#   |                     from the literature
#   |
#   |-- 2_ENV_VARIABLES/  Raster files (.tif) with the environmental variables.
#   |                     The rasters are stacked in alphabetical order of file
#   |                     name and the layer names set in the LOAD ENVIRONMENTAL
#   |                     VARIABLES section follow that order, so the folder must
#   |                     hold exactly the 19 predictor rasters and no other
#   |                     files matching ".tif" (e.g. .tif.aux.xml)
#   |
#   |-- 3_STUDY_AREA/     Shapefile of the study area (.shp, with its
#   |                     .dbf, .shx and .prj files)
#   |
#   `-- 4_RESULTS/        Output folder for the modelling results (correlation
#                         matrix, VIF tables, selected variables, figures, ...).
#                         Create it before running the script: the script does
#                         not create it.
#
# Besides 4_RESULTS, biomod2 creates a folder named after the species
# ("Alouatta.belzebul") inside the project folder, with the fitted models and
# projections, and the ensemble tables (EM*.csv) are saved in the project folder.


# ==============================================================================
# PACKAGES
# ==============================================================================
# Package versions: CRAN snapshot of 2024-06-01 (e.g. biomod2 4.2-5).

# ---- Version check and installation ------------------------------------------
# Every package (and each of its dependencies) is compared with the version that
# was current on 2024-06-01, and installed from a dated copy of CRAN (Posit Public
# Package Manager) only when it is missing or in a different version. If
# everything already matches, nothing is installed.
# - Needs an internet connection to read the snapshot index; without it, the
#   check is skipped with a warning.
# - If packages have to be installed, run this script in a fresh R session
#   (Session > Restart R), before any package is loaded.
# - Packages that ship with R (base and recommended, e.g. stats, MASS, Matrix)
#   are not touched.

snapshot_date <- "2024-06-01"
snapshot_repo <- paste0("https://packagemanager.posit.co/cran/", snapshot_date)

# Packages loaded by this script, with the version expected in the snapshot
used_versions <- c(
  sf = "1.0-16", ggplot2 = "3.5.1", raster = "3.6-26", ggcorrplot = "0.1.4.1",
  viridis = "0.6.5", dplyr = "1.1.4", tidyr = "1.3.1", usdm = "2.1-7",
  terra = "1.7-78", dismo = "1.3-14", biomod2 = "4.2-5", wesanderson = "0.3.7",
  ggthemes = "5.1.0"
)

# Packages that biomod2 calls internally for the algorithms used below
# (RF, ANN, MARS, FDA, CTA, XGBOOST). They are only installed, not loaded.
backend_packages <- c("randomForest", "nnet", "earth", "mda", "rpart", "xgboost")

#' Install the snapshot version of packages and of all their dependencies
#'
#' @param pkgs     Character vector of package names.
#' @param expected Named character vector with the versions the script was
#'                 written for. A different version gives a warning, except for
#'                 biomod2, which stops the script.
#' @param repo     URL of the dated CRAN snapshot.
install_snapshot_packages <- function(pkgs, expected, repo = snapshot_repo) {
  
  # Index of the snapshot (package names and versions available on that date)
  avail <- available.packages(repos = repo, type = "source")
  if (nrow(avail) == 0) {
    warning("Could not read the snapshot index (no internet connection?). ",
            "Package version check skipped.", call. = FALSE)
    return(invisible(FALSE))
  }
  
  not_in_snapshot <- setdiff(pkgs, rownames(avail))
  if (length(not_in_snapshot) > 0) {
    stop("Not found in the snapshot (", snapshot_date, "): ",
         paste(not_in_snapshot, collapse = ", "), call. = FALSE)
  }
  
  # Requested packages + all recursive dependencies, minus what ships with R
  shipped_with_r <- c("R", rownames(installed.packages(priority = "high")))
  deps <- tools::package_dependencies(pkgs, db = avail, recursive = TRUE,
                                      which = c("Depends", "Imports", "LinkingTo"))
  all_pkgs <- setdiff(unique(c(pkgs, unlist(deps))), shipped_with_r)
  
  # TRUE when the package is installed in exactly the snapshot version.
  # package_version() treats "4.2-5" and "4.2.5" as the same version.
  is_current <- function(pkg) {
    if (!nzchar(system.file(package = pkg))) return(FALSE)
    packageVersion(pkg) == package_version(avail[pkg, "Version"])
  }
  to_install <- all_pkgs[!vapply(all_pkgs, is_current, logical(1))]
  
  if (length(to_install) == 0) {
    message("[install] All packages already match the ", snapshot_date, " snapshot.")
  } else {
    # Files of a package that is loaded in this session cannot be replaced safely
    # (Windows locks them), so ask for a fresh session first.
    loaded <- intersect(to_install, loadedNamespaces())
    if (length(loaded) > 0) {
      stop("These packages are loaded in this session and need another version: ",
           paste(loaded, collapse = ", "),
           ".\nRestart R (Session > Restart R) and run the script again.",
           call. = FALSE)
    }
    message("[install] Installing ", length(to_install), " package(s) from the ",
            snapshot_date, " snapshot:\n  ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = repo)
  }
  
  # Check that everything ended up in the snapshot version
  failed <- all_pkgs[!vapply(all_pkgs, is_current, logical(1))]
  if (length(failed) > 0) {
    stop("Could not install the snapshot version of: ",
         paste(failed, collapse = ", "), call. = FALSE)
  }
  
  # Compare with the versions this script was written for
  for (pkg in names(expected)) {
    if (packageVersion(pkg) != package_version(expected[[pkg]])) {
      msg <- paste0(pkg, " is version ", packageVersion(pkg), " but ",
                    expected[[pkg]], " was expected.")
      if (pkg == "biomod2") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
    }
  }
  invisible(TRUE)
}

install_snapshot_packages(c(names(used_versions), backend_packages), used_versions)

# ---- Load packages -----------------------------------------------------------
# Loading order matters because of function masking (e.g. dplyr::select vs raster::select)

library(sf) # Used version: 1.0-16
library(ggplot2) # Used version: 3.5.1
library(raster) # Used version: 3.6-26
library(ggcorrplot) # Used version: 0.1.4.1
library(viridis) # Used version: 0.6.5
library(dplyr) # Used version: 1.1.4
library(tidyr) # Used version: 1.3.1
library(usdm) # Used version: 2.1-7
library(terra) # Used version: 1.7-78
library(dismo) # Used version: 1.3-14
library(biomod2) # Used version: 4.2-5
library(wesanderson) # Used version: 0.3.7
library(ggthemes) # Used version: 5.1.0


# ==============================================================================
# SETTING SEED
# ==============================================================================
set.seed(2222)  # reproducible random sampling (background points, etc.)


# ==============================================================================
# WORKING DIRECTORY
# ==============================================================================
# Edit this path to match your machine. All relative paths below ("1_OCCURRENCES",
# "2_ENV_VARIABLES", etc.) are resolved from here.
setwd('C:/path/to/your/working/directory')


# ==============================================================================
# LOAD OCCURRENCE DATA
# ==============================================================================
# Literature records (shapefile)
occ_lit <- read_sf('1_OCCURRENCES/pontos_abelzebul_lit.shp')

# Field sampling records (CSV converted to an sf object, WGS84).
# Not used in the modelling steps below yet: reserved for model validation (A2).
occ_samp <- read.csv("1_OCCURRENCES/centroides_detec_alouatta_belzebul.csv") %>% 
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 'EPSG:4326')


# ==============================================================================
# LOAD STUDY AREA
# ==============================================================================
study_area <- read_sf('3_STUDY_AREA/buffer_cep.shp')


# ==============================================================================
# LOAD ENVIRONMENTAL VARIABLES
# ==============================================================================
# Stack all GeoTIFFs in the folder into a single multi-layer raster
env_variables <- stack(list.files('2_ENV_VARIABLES', pattern = ".tif", full.names = TRUE))

# Restrict the rasters to the study area
env_variables <- mask(env_variables, study_area)

# Layer names.
names(env_variables) <- c("canopyHeight", "treeCanopyCover", "slope", "distAreasProt", "distRoads", "elevation", "gHM",
                          "nFrag", "percAgropastoral", "percOtherForests", "annualPrecip", "precDryMonth", "precSeason",
                          "precWetMonth", "tempMaxWarmMonth", "annualTemp", "tempMinColdMonth", "tempAnnualRange",
                          "tempSeason")


# ==============================================================================
# PREDICTOR SELECTION: CORRELATION AND COLLINEARITY
# ==============================================================================

# ---- Step 1: Pairwise correlation (Pearson) ----------------------------------
# Correlation matrix across all raster layers
corMatrix <- as.data.frame(layerStats(env_variables, 'pearson', na.rm = T)$`pearson correlation coefficient`)

# Correlogram (variables ordered by hierarchical clustering)
ggcorrplot(corMatrix, hc.order = TRUE, lab = TRUE, ggtheme = ggplot2::theme_light,
           show.legend = TRUE, legend.title = "Pearson's R", show.diag = FALSE,
           colors = inferno(3),
           outline.color = "white",
           lab_size = 3, pch = .1, pch.col = "black", pch.cex = .1,
           tl.cex = 5, digits = 2) +
  theme(axis.text.x = element_text(size = 11, angle = 90, hjust = 0.95, vjust = 0.5), 
        axis.text.y = element_text(size = 11, hjust = 0.95, vjust = 0.5))

# ggsave() saves the last plot that was drawn
ggsave("4_RESULTS/correlationMatrix.pdf", scale = 1.25, dpi = 600)

write.csv(corMatrix, "4_RESULTS/correlationMatrix.csv")

# List of highly correlated pairs (|r| > 0.7). The upper bound (0.999999) removes
# the diagonal (a variable correlated with itself).
corVariables <- layerStats(env_variables, 'pearson', na.rm = T)$`pearson correlation coefficient` %>%
  as.data.frame() %>%
  mutate(var1 = rownames(.)) %>%
  gather(var2, value, -var1) %>%
  arrange(desc(value)) %>%
  filter(abs(value) > 0.7 & abs(value) < 0.999999)

# Number of times each variable appears in a highly correlated pair
# (helps decide which variables to drop first)
corVariablesNumber <- as.data.frame(table(corVariables$var1)[order(table(corVariables$var1), decreasing = TRUE)])

# Variables removed after inspecting the correlation results
toRemove <- c("nFrag", "precWetMonth", "precDryMonth", "tempAnnualRange", "annualTemp", "tempMinColdMonth")

env_variables <- dropLayer(env_variables, toRemove)

plot(env_variables)

# ---- Step 2: Multicollinearity (VIF) -----------------------------------------
# VIF of the remaining variables (usdm::vif works on terra rasters)
vifCollinearity <- vif(rast(env_variables))
write.csv(vifCollinearity, '4_RESULTS/vifCollinearity_raw.csv')

# Rule: among variables with VIF > 4, remove one at a time the one with the
# lowest ecological relevance, then recompute VIF.
toRemove2 <- c("elevation")

vifCollinearity <- vif(rast(dropLayer(env_variables, toRemove2)))
write.csv(vifCollinearity, '4_RESULTS/vifCollinearity_final.csv')

env_variables <- dropLayer(env_variables, toRemove2)
plot(env_variables)

# Save the final set of predictors
terra::writeRaster(rast(env_variables), "4_RESULTS/selectedVariables.tif", overwrite = TRUE)


# ==============================================================================
# MODELLING
# ==============================================================================

# ---- Occurrence data ---------------------------------------------------------
# Keep only one occurrence per raster pixel (reduces spatial pseudoreplication)
s <- gridSample(data.frame(x = occ_lit$Longitude, y = occ_lit$Latitude), env_variables[[1]], n = 1) 
occ_lit <- occ_lit[row.names(s), ]

# Presence column (1 = presence)
occ_lit$Presenca <- 1

# Rename columns
colnames(occ_lit) <- c("Species", "Method", "Locality", "Municipality", "State", "Start_date", "End_date", "Longitude", "Latitude", "Imprecision", "Datum",
                       "Ref_ID", "Ref_type", "Reference", "Source_database", "field_16", "geometry", "Presence")

# Replace spaces in species names with underscores (e.g. "Alouatta_belzebul")
occ_lit$Species <- gsub(" ", "_", occ_lit$Species)

# ---- Background points -------------------------------------------------------
# Draw 10,000 random points inside the study area
back_data <- st_sample(
  study_area,
  10000,
  type = "random",
  exact = TRUE,
  warn_if_not_integer = TRUE,
  by_polygon = FALSE,
  progress = TRUE) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  st_as_sf(coords = c("X", "Y"), remove = FALSE)

# Give the background table the same columns as the presence table
# (all fields NA except species, coordinates and geometry)
back_data <- data.frame("Species" = unique(occ_lit$Species), "Method" = NA, "Locality" = NA, "Municipality" = NA, "State" = NA, "Start_date" = NA, "End_date" = NA,
                        "Longitude" = back_data$X, "Latitude" = back_data$Y, "Imprecision" = NA, "Datum" = NA, "Ref_ID" = NA, "Ref_type" = NA,
                        "Reference" = NA, "Source_database" = NA, "field_16" = NA, "geometry" = back_data$geometry, "Presence" = NA)

# Combine presences and background points
belzebul <- rbind(as.data.frame(occ_lit), back_data)


# ---- Format data for biomod2 -------------------------------------------------
myBiomodData <- BIOMOD_FormatingData(resp.var = belzebul$Presence,
                                     expl.var = stack(env_variables),
                                     resp.xy = belzebul[, c("Longitude", "Latitude")],
                                     resp.name = unique(belzebul$Species))

# ---- Environmental values at presence and pseudo-absence points --------------
dataTable <- as.data.frame(myBiomodData@data.env.var)
dataTable <- cbind(myBiomodData@coord, dataTable)
dataTable$Presence <- myBiomodData@data.species
dataTable$Presence <- ifelse(is.na(dataTable$Presence), "Pseudo-absence", "Presence")
dataTable$Species <- "Alouatta belzebul"

# Move the species column to the front.
# Columns are referenced by position: 2 coordinates + 12 predictors + Presence
# + Species = 16. Update the indices if the number of predictors changes.
dataTable <- dataTable[, c(16, 1:15)]

write.csv(dataTable, "4_RESULTS/myBiomodData.csv")

# ---- Distribution of each predictor: presences vs pseudo-absences ------------
# Columns 4:15 are the 12 predictors; reshape to long format for plotting
dataTable2 <- dataTable[, 4:15] %>%
  gather(key = "text", value = "value") %>%
  mutate(text = gsub("\\.", " ", text)) %>%
  mutate(value = as.numeric(value))
dataTable2$presence <- dataTable$Presence

gg1 <- ggplot(dataTable2, aes(x = value, fill = presence)) +
  geom_density(colour = "black", alpha = .7, linewidth = .1) +
  scale_fill_manual(values = wes_palette(n = 3, name = 'Darjeeling1'),
                    breaks = c("Pseudo-absence", "Presence")) +
  theme_light() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text = element_text(size = 5)) +
  xlab("") + ylab("") +
  facet_wrap(~ text, scales = 'free')

gg1

ggsave("4_RESULTS//gg1_distribuicao_pres_aus.png", gg1, dpi = 600)


# ==============================================================================
# SINGLE MODELS
# ==============================================================================
# 8 algorithms x 100 random 80/20 calibration/validation splits.
# var.import = 0 skips variable importance here (it is computed for the ensemble).
# nb.cpu = number of cores; change it to match your machine.
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    models = c("GLM", "RF", "ANN", "MARS", "FDA", "CTA", "SRE", "XGBOOST"),
                                    CV.nb.rep = 100,
                                    data.split.perc = 80,
                                    var.import = 0,
                                    metric.eval = c('TSS'),
                                    CV.do.full.models = FALSE,
                                    prevalence = NULL,
                                    seed.val = 2222,
                                    nb.cpu = 18)

# ---- Evaluation plots --------------------------------------------------------
# gg2: evaluation scores by algorithm (dotted line = TSS threshold of 0.8)
gg2 <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
gg2_1 <- gg2$plot + theme_light() +
  scale_fill_brewer(palette = 'Accent') +
  geom_hline(data = data.frame(yint = 0.8, Eval.metric = "TSS"), aes(yintercept = yint), linetype = "dotted") +
  xlab("") + ylab("") +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

gg2_1

# gg3: evaluation scores by cross-validation run
gg3 <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('run', 'run'))
gg3 <- gg3$plot + theme_light() +
  scale_color_viridis(discrete = T, option = 'plasma') +
  geom_hline(data = data.frame(yint = 0.8, Eval.metric = "TSS"), aes(yintercept = yint), linetype = "dotted") +
  xlab("") + ylab("") +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
gg3


ggsave("4_RESULTS//gg2_1_aval_alg.png", gg2_1, dpi = 600, width = 4199, height = 3199, units = "px")
ggsave("4_RESULTS//gg3_aval_runs.png", gg3, dpi = 600, width = 19000, height = 8000, units = "px", limitsize = F)


# ==============================================================================
# ENSEMBLE MODELS
# ==============================================================================
# Reload the saved single-model output from disk (useful when the single models
# were fitted in a previous session). The path is relative to the working
# directory used when BIOMOD_Modeling() was run.
# myBiomodModelOut <- assign("myBiomodModelOut", get(load("Alouatta.belzebul/.BIOMOD_DATA/AllModels/Alouatta.belzebul.AllModels.models.out")))


# Weighted-mean ensemble (EMwmean) of all single models with TSS >= 0.8.
# Weights are proportional to each model's TSS.
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.8),
                                      var.import = 3,
                                      metric.eval = c('TSS'),
                                      em.algo = 'EMwmean',
                                      EMwmean.decay = 'proportional',
                                      seed.val = 2222)

# ---- Evaluation scores of the ensemble ---------------------------------------
EMevaluations <- get_evaluations(myBiomodEM)
write.csv(EMevaluations, "EMevaluations.csv")

EMevaluations$algo <- "EM"

# ---- Variable importance of the ensemble -------------------------------------
# Per run (3 permutations)
EMvarImportance <- get_variables_importance(myBiomodEM)
write.csv(EMvarImportance, "EMvarImportance_run.csv")

# Mean importance per variable, sorted from most to least important
EMvarImportance <- EMvarImportance %>%
  group_by(expl.var) %>%
  summarize(m = mean(var.imp)) %>%
  arrange(desc(m), expl.var)
write.csv(EMvarImportance, "EMvarImportance_mean.csv")

# ---- Plots with the ensemble -------------------------------------------------
# gg2_2: same as gg2, with the ensemble ("EM") added as one more "algorithm"
gg2_2 <- gg2

# Ensemble scores in the same format as the single-model evaluation table
EMevaluations2 <- data.frame('full.name' = EMevaluations$full.name, 
                             'PA' = 'allData',
                             'run' = 'RUN1',
                             'algo' = EMevaluations$algo,
                             'metric.eval' = EMevaluations$metric.eval,
                             'cutoff' = EMevaluations$cutoff,
                             'sensitivity' = as.numeric(EMevaluations$sensitivity),
                             'specificity' = as.numeric(EMevaluations$specificity),
                             'calibration' = as.numeric(EMevaluations$calibration),
                             'validation' = NA,
                             'evaluation' = NA)

# Append the ensemble to the single-model plot data and redraw
gg2_2$plot$data <- rbind(gg2_2$plot$data, EMevaluations2)
gg2_2 <- gg2_2$plot + theme_light() + 
  geom_hline(data = data.frame(yint = 0.8, Eval.metric = "TSS"), aes(yintercept = yint), linetype = "dotted") +
  scale_fill_brewer(palette = 'Accent')  +
  xlab("") + ylab("") +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
gg2    # single models only
gg2_2  # single models + ensemble

# gg4: mean variable importance in the ensemble
gg4 <- ggplot() +
  geom_bar(data = EMvarImportance, aes(x = reorder(expl.var, m), y = m, fill = reorder(expl.var, m)), stat = "identity") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.7, direction = 1, option = 'inferno') +
  coord_flip() +
  theme_light() +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        plot.title = element_text(size = 11),
        plot.margin = unit(c(0.1, 0.5, 0, -0.3), "cm")) +
  xlab("") + ylab("") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

gg4

# gg5: response curves of the ensemble (other predictors fixed at their median)
gg5 <- bm_PlotResponseCurves(bm.out = myBiomodEM,
                             models.chosen = 'all',
                             fixed.var = 'median')$plot

# Keep only the curves of the final ensemble (weighted mean by TSS)
gg5$data <- subset(gg5$data, pred.name == "Alouatta.belzebul_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo")

gg5 <- gg5 + theme_light() + theme(legend.position = "none") +
  geom_line(aes(color = expl.name)) +
  scale_color_viridis(discrete = TRUE, alpha = 1, direction = -1, option = 'inferno') +
  theme_light() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.3), "cm")) +
  facet_wrap(~expl.name, scales = 'free', ncol = 3) +
  labs(x = '', y = 'Occurrence probability') + ggtitle("")

gg5

# ---- Save figures ------------------------------------------------------------
# Each figure is saved twice: in the results folder (PNG) and in the manuscript
# folder (PDF)
ggsave("4_RESULTS//gg2_2.png", gg2_2, dpi = 600, width = 4199, height = 3199, units = "px")

ggsave("4_RESULTS//gg4.png", gg4, dpi = 600, width = 4199, height = 4199, units = "px")

ggsave("4_RESULTS//gg5.png", gg5, dpi = 600, width = 4199, height = 4199, units = 'px')

# Save the response-curve data
write.csv(gg5$data, "EMresponseCurves.csv")


# ==============================================================================
# PROJECTION ONTO CURRENT CONDITIONS
# ==============================================================================
# To reload a saved ensemble model from a previous session instead:
# myBiomodEM <- assign("myBiomodEM", get(load(paste0(getwd(), "/", sp, "/", sp, ".AllModels.ensemble.models.out"))))

# Project the ensemble over the study area. Also produces binary maps (using the
# TSS-optimised threshold), filtered maps, and a clamping mask showing where the
# environmental values fall outside the range used for calibration.
myBiomodProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                           proj.name = 'CurrentEM',
                                           new.env = env_variables,
                                           models.chosen = 'all',
                                           metric.binary = c('TSS'),
                                           metric.filter = c('TSS'),
                                           build.clamping.mask = TRUE,
                                           output.format = '.img',
                                           seed.val = 2222)
