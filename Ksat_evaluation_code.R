############################################################
# Validation of Ksat Map (0–20 cm)
############################################################

library(dplyr)
library(terra)
library(tidyr)
library(ggplot2)

download.file(
  "https://R-Forge.R-project.org/src/contrib/ithir_1.0.tar.gz",
  "ithir.tar.gz"
)

untar("ithir.tar.gz")

list.files(
  "ithir/R",
  pattern = "ea_spline",
  full.names = TRUE
)

source("ithir/R/ea_spline.R")

exists("ea_spline")

############################################################
# Load EUHYDI datasets
############################################################

general_data <- read.csv(
  "C:/Users/surya/Downloads/EUHYDI_public_v1.1.0_csv/GENERAL.csv"
)

general_data <- general_data[, 1:7]

basic_data <- read.csv(
  "C:/Users/surya/Downloads/EUHYDI_public_v1.1.0_csv/Basic.csv"
)

conductivity_data <- read.csv(
  "C:/Users/surya/Downloads/EUHYDI_public_v1.1.0_csv/Cond.csv"
)

############################################################
# Select saturated hydraulic conductivity (HEAD = 0)
############################################################

ksat_measurements <- conductivity_data %>%
  filter(VALUE == 0)

ksat_measurements <- ksat_measurements[, 1:6]

############################################################
# Average duplicate measurements
############################################################

ksat_measurements_mean <- aggregate(
  . ~ SAMPLE_ID,
  data = ksat_measurements,
  mean
)

############################################################
# Merge profile information
############################################################

ksat_profile_data <- merge(
  general_data,
  ksat_measurements_mean,
  by = "PROFILE_ID"
)

ksat_profile_basic <- merge(
  ksat_profile_data,
  basic_data,
  by = "SAMPLE_ID"
)

############################################################
# Keep required columns
############################################################

ksat_profile_coordinates <- ksat_profile_basic[, c(
  1, 2,
  6, 7,
  10, 11,
  15, 16
)]

ksat_profile_coordinates <- ksat_profile_coordinates %>%
  filter(!is.na(X_WGS84))

############################################################
# Prepare spline input
############################################################

ksat_spline_input <- ksat_profile_coordinates[, c(
  "PROFILE_ID.x",
  "SAMPLE_DEP_TOP",
  "SAMPLE_DEP_BOT",
  "COND"
)]

colnames(ksat_spline_input) <- c(
  "PROFILE_ID",
  "TOP",
  "BOT",
  "COND"
)

############################################################
# Remove invalid horizons
############################################################

ksat_spline_input <- ksat_spline_input %>%
  filter(
    TOP >= 0,
    BOT > TOP,
    BOT > 0
  ) %>%
  filter(
    complete.cases(.)
  )

############################################################
# Average duplicate horizons
############################################################

ksat_spline_input <- ksat_spline_input %>%
  group_by(PROFILE_ID, TOP, BOT) %>%
  summarise(
    COND = mean(COND, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
# Keep profiles having at least two horizons
############################################################

ksat_spline_input <- ksat_spline_input %>%
  group_by(PROFILE_ID) %>%
  filter(n() >= 2) %>%
  ungroup()

############################################################
# Harmonize to standard depths
############################################################

ksat_spline <- ea_spline(
  obj = ksat_spline_input,
  var.name = "COND",
  d = c(0, 20, 30, 60, 100, 200),
  lam = 0.1,
  vlow = 0,
  show.progress = FALSE
)

############################################################
# Extract 0–20 cm Ksat
############################################################

ksat_harmonised <- ksat_spline$harmonised

ksat_depth020 <- ksat_harmonised[, c(
  "id",
  "0-20 cm"
)]

colnames(ksat_depth020)[1] <- "PROFILE_ID.x"

############################################################
# Extract profile coordinates
############################################################

ksat_coordinates <- ksat_profile_coordinates[, c(
  "PROFILE_ID.x",
  "X_WGS84",
  "Y_WGS84"
)]

############################################################
# Merge coordinates with harmonized values
############################################################

ksat_profiles <- merge(
  ksat_depth020,
  ksat_coordinates,
  by = "PROFILE_ID.x"
)

############################################################
# Keep one point per profile
############################################################

ksat_profiles_unique <- ksat_profiles %>%
  distinct(PROFILE_ID.x, .keep_all = TRUE)


############################################################
# Load Ksat prediction raster
############################################################

ksat_raster <- rast(
  "C:/Users/surya/Downloads/Panos_SHP_data/Resample_1km/Resample_Ks_PTF02_modified_6_03_tif_non_log.tif"
)

############################################################
# Create spatial points
############################################################

ksat_points <- vect(
  ksat_profiles_unique,
  geom = c("X_WGS84", "Y_WGS84"),
  crs = "EPSG:4326"
)

############################################################
# Extract predicted Ksat values
############################################################

ksat_prediction <- extract(
  ksat_raster,
  ksat_points
)

ksat_profiles_unique$Predicted_Ksat <- ksat_prediction[,2]

############################################################
# Remove missing values
############################################################

ksat_validation <- ksat_profiles_unique %>%
  filter(
    complete.cases(.)
  )

############################################################
# Log10 transformation
############################################################

ksat_validation$Measured_log10 <-
  log10(ksat_validation$`0-20 cm`)

ksat_validation$Predicted_log10 <-
  log10(ksat_validation$Predicted_Ksat)

############################################################
# Remove extremely small measured values
############################################################

ksat_validation <- ksat_validation %>%
  filter(Measured_log10 > -6)

############################################################
# Model performance
############################################################

MAE <- mean(
  abs(
    ksat_validation$Measured_log10 -
      ksat_validation$Predicted_log10
  )
)

RMSE <- sqrt(
  mean(
    (
      ksat_validation$Measured_log10 -
        ksat_validation$Predicted_log10
    )^2
  )
)

Bias <- mean(
  ksat_validation$Predicted_log10 -
    ksat_validation$Measured_log10
)

############################################################
# Print statistics
############################################################

cat("MAE  =", round(MAE,3), "\n")
cat("RMSE =", round(RMSE,3), "\n")
cat("Bias =", round(Bias,3), "\n")


############################################################
# Load EU-SoilHydroGrids raster
############################################################

ksat_euhydro_raster <- rast(
  "X:/Users/Surya/EU_SoilHydroGrids_1km/Ks_0_20cm_mask_log10.tif"
)

############################################################
# Create spatial points
############################################################

ksat_points <- vect(
  ksat_validation,
  geom = c("X_WGS84","Y_WGS84"),
  crs = "EPSG:4326"
)

############################################################
# Extract EU-SoilHydroGrids predictions
############################################################

ksat_euhydro_prediction <- extract(
  ksat_euhydro_raster,
  ksat_points
)

ksat_validation$EUHydroGrid_log10 <-
  ksat_euhydro_prediction[,2]

ksat_validation <- ksat_validation %>%
  filter(
    complete.cases(.)
  )

############################################################
# Model performance
############################################################


MAE <- mean(
  abs(
    ksat_validation$Measured_log10 -
      ksat_validation$EUHydroGrid_log10
  )
)

RMSE <- sqrt(
  mean(
    (
      ksat_validation$Measured_log10 -
        ksat_validation$EUHydroGrid_log10
    )^2
  )
)

Bias <- mean(
  ksat_validation$EUHydroGrid_log10 -
    ksat_validation$Measured_log10
)

############################################################
# Print statistics
############################################################

cat("MAE  =", round(MAE,3), "\n")
cat("RMSE =", round(RMSE,3), "\n")
cat("Bias =", round(Bias,3), "\n")

############################################################
# Density plots
############################################################

density_data <- bind_rows(
  
  data.frame(
    Value = ksat_validation$Measured_log10,
    Type = "Measured",
    Map = "EUPTFv2–LUCAS map"
  ),
  
  data.frame(
    Value = ksat_validation$Predicted_log10,
    Type = "Predicted",
    Map = "EUPTFv2–LUCAS map"
  ),
  
  data.frame(
    Value = ksat_validation$Measured_log10,
    Type = "Measured",
    Map = "EU-SoilHydroGrids map"
  ),
  
  data.frame(
    Value = ksat_validation$EUHydroGrid_log10,
    Type = "Predicted",
    Map = "EU-SoilHydroGrids map"
  )
  
)

ggplot(
  density_data,
  aes(
    x = Value,
    fill = Type,
    colour = Type
  )
) +
  geom_density(
    alpha = 0.3,
    linewidth = 1
  ) +
  facet_wrap(~Map) +
  theme_bw(base_size = 14) +
  labs(
    x = expression(log[10](Ksat)),
    y = "Density"
  )
