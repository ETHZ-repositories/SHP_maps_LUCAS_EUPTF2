
############################################################
# Validation of Permanent Wilting Point (θ1500)
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

retention_data <- read.csv(
  "C:/Users/surya/Downloads/EUHYDI_public_v1.1.0_csv/RET.csv"
)

############################################################
# Select permanent wilting point measurements
# (Closest measurement to 15000 cm pressure head)
############################################################

wp_measurements <- retention_data %>%
  filter(HEAD >= 15000 & HEAD <= 16000)

wp_measurements <- wp_measurements %>%
  group_by(SAMPLE_ID) %>%
  slice_min(
    abs(HEAD - 15000),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

wp_measurements <- wp_measurements[, 1:4]

############################################################
# Average duplicate measurements
############################################################

wp_measurements_mean <- aggregate(
  . ~ SAMPLE_ID,
  data = wp_measurements,
  mean
)

############################################################
# Merge profile information
############################################################

wp_profile_data <- merge(
  general_data,
  wp_measurements_mean,
  by = "PROFILE_ID"
)

wp_profile_basic <- merge(
  wp_profile_data,
  basic_data,
  by = "SAMPLE_ID"
)

############################################################
# Keep required columns
############################################################

wp_profile_coordinates <- wp_profile_basic[, c(
  1, 2,
  6, 7,
  9, 10,
  13, 14
)]

wp_profile_coordinates <- wp_profile_coordinates %>%
  filter(!is.na(X_WGS84))

############################################################
# Prepare spline input
############################################################

wp_spline_input <- wp_profile_coordinates[, c(
  "PROFILE_ID.x",
  "SAMPLE_DEP_TOP",
  "SAMPLE_DEP_BOT",
  "THETA"
)]

colnames(wp_spline_input) <- c(
  "PROFILE_ID",
  "TOP",
  "BOT",
  "THETA"
)

############################################################
# Remove invalid horizons
############################################################

wp_spline_input <- wp_spline_input %>%
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

wp_spline_input <- wp_spline_input %>%
  group_by(PROFILE_ID, TOP, BOT) %>%
  summarise(
    THETA = mean(THETA, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
# Keep profiles with at least two horizons
############################################################

wp_spline_input <- wp_spline_input %>%
  group_by(PROFILE_ID) %>%
  filter(n() >= 2) %>%
  ungroup()

############################################################
# Harmonize to standard depths
############################################################

wp_spline <- ea_spline(
  obj = wp_spline_input,
  var.name = "THETA",
  d = c(0, 20, 30, 60, 100, 200),
  lam = 0.1,
  vlow = 0,
  show.progress = FALSE
)

############################################################
# Extract 0–20 cm permanent wilting point
############################################################

wp_harmonised <- wp_spline$harmonised

wp_depth020 <- wp_harmonised[, c(
  "id",
  "0-20 cm"
)]

colnames(wp_depth020)[1] <- "PROFILE_ID.x"

############################################################
# Extract profile coordinates
############################################################

wp_coordinates <- wp_profile_coordinates[, c(
  "PROFILE_ID.x",
  "X_WGS84",
  "Y_WGS84"
)]

############################################################
# Merge coordinates with harmonized values
############################################################

wp_profiles <- merge(
  wp_depth020,
  wp_coordinates,
  by = "PROFILE_ID.x"
)

############################################################
# Keep one point per profile
############################################################

wp_profiles_unique <- wp_profiles %>%
  distinct(PROFILE_ID.x, .keep_all = TRUE)
######################################################


############################################################
# Load Wilting Point prediction raster
############################################################

wp_raster <- rast(
  "C:/Users/surya/Downloads/Panos_SHP_data/Resample_1km/Resample_WP_150kpa_PTF07_modified_6_03_tif.tif"
)

############################################################
# Create spatial points
############################################################

wp_points <- vect(
  wp_profiles_unique,
  geom = c("X_WGS84", "Y_WGS84"),
  crs = "EPSG:4326"
)

############################################################
# Extract predicted Wilting Point values
############################################################

wp_prediction <- extract(
  wp_raster,
  wp_points
)

wp_profiles_unique$Predicted_WP <- wp_prediction[, 2]

############################################################
# Remove missing values
############################################################

wp_validation <- wp_profiles_unique %>%
  filter(
    complete.cases(.)
  )

############################################################
# Model performance
############################################################

MAE <- mean(
  abs(
    wp_validation$`0-20 cm` -
      wp_validation$Predicted_WP
  )
)

RMSE <- sqrt(
  mean(
    (
      wp_validation$`0-20 cm` -
        wp_validation$Predicted_WP
    )^2
  )
)

Bias <- mean(
  wp_validation$Predicted_WP -
    wp_validation$`0-20 cm`
)

############################################################
# Print statistics
############################################################

cat("MAE  =", round(MAE, 3), "\n")
cat("RMSE =", round(RMSE, 3), "\n")
cat("Bias =", round(Bias, 3), "\n")

###################################

############################################################
# Load EU-SoilHydroGrids raster
############################################################

wp_euhydro_raster <- rast(
  "X:/Users/Surya/EU_SoilHydroGrids_1km/WP_0_20cm_mask.tif"
)

############################################################
# Create spatial points
############################################################

wp_points <- vect(
  wp_validation,
  geom = c("X_WGS84", "Y_WGS84"),
  crs = "EPSG:4326"
)

############################################################
# Extract EU-SoilHydroGrids predictions
############################################################

wp_euhydro_prediction <- extract(
  wp_euhydro_raster,
  wp_points
)

wp_validation$EUHydroGrid_WP <-
  wp_euhydro_prediction[, 2] / 100

############################################################
# Remove missing values
############################################################

wp_validation <- wp_validation %>%
  filter(
    complete.cases(.)
  )

############################################################
# Model performance
############################################################

MAE <- mean(
  abs(
    wp_validation$`0-20 cm` -
      wp_validation$EUHydroGrid_WP
  )
)

RMSE <- sqrt(
  mean(
    (
      wp_validation$`0-20 cm` -
        wp_validation$EUHydroGrid_WP
    )^2
  )
)

Bias <- mean(
  wp_validation$EUHydroGrid_WP -
    wp_validation$`0-20 cm`
)

############################################################
# Print statistics
############################################################

cat("MAE  =", round(MAE, 3), "\n")
cat("RMSE =", round(RMSE, 3), "\n")
cat("Bias =", round(Bias, 3), "\n")

###############################################################
plot_data <- bind_rows(
  data.frame(
    Value = wp_validation$`0-20 cm`,
    Type = "Measured",
    Map = "EUPTFv2–LUCAS map"
  ),
  data.frame(
    Value = wp_validation$Predicted_WP,
    Type = "Predicted",
    Map = "EUPTFv2–LUCAS map"
  ),
  data.frame(
    Value = wp_validation$`0-20 cm`,
    Type = "Measured",
    Map = "EU-SoilHydroGrids map"
  ),
  data.frame(
    Value = wp_validation$EUHydroGrid_WP,
    Type = "Predicted",
    Map = "EU-SoilHydroGrids map"
  )
)

ggplot(plot_data,
       aes(x = Value,
           fill = Type,
           colour = Type)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  facet_wrap(~Map) +
  theme_bw(base_size = 14) +
  labs(
    x = expression(WP~"["*cm^3/cm^3*"]"),
    y = "Density"
  )
