
############################################################
# Validation of Field Capacity (θ33)
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
# Select field capacity measurements
# (Closest measurement to 330 cm pressure head)
############################################################

fc_measurements <- retention_data %>%
  filter(HEAD >= 300 & HEAD <= 350)

fc_measurements <- fc_measurements %>%
  group_by(SAMPLE_ID) %>%
  slice_min(
    abs(HEAD - 330),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

fc_measurements <- fc_measurements[, 1:4]

############################################################
# Average duplicate measurements
############################################################

fc_measurements_mean <- aggregate(
  . ~ SAMPLE_ID,
  data = fc_measurements,
  mean
)

############################################################
# Merge profile information
############################################################

fc_profile_data <- merge(
  general_data,
  fc_measurements_mean,
  by = "PROFILE_ID"
)

fc_profile_basic <- merge(
  fc_profile_data,
  basic_data,
  by = "SAMPLE_ID"
)

############################################################
# Keep required columns
############################################################

fc_profile_coordinates <- fc_profile_basic[, c(
  1, 2,
  6, 7,
  9, 10,
  13, 14
)]

fc_profile_coordinates <- fc_profile_coordinates %>%
  filter(!is.na(X_WGS84))

############################################################
# Prepare spline input
############################################################

fc_spline_input <- fc_profile_coordinates[, c(
  "PROFILE_ID.x",
  "SAMPLE_DEP_TOP",
  "SAMPLE_DEP_BOT",
  "THETA"
)]

colnames(fc_spline_input) <- c(
  "PROFILE_ID",
  "TOP",
  "BOT",
  "THETA"
)

############################################################
# Remove invalid horizons
############################################################

fc_spline_input <- fc_spline_input %>%
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

fc_spline_input <- fc_spline_input %>%
  group_by(PROFILE_ID, TOP, BOT) %>%
  summarise(
    THETA = mean(THETA, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
# Keep profiles with at least two horizons
############################################################

fc_spline_input <- fc_spline_input %>%
  group_by(PROFILE_ID) %>%
  filter(n() >= 2) %>%
  ungroup()

############################################################
# Harmonize to standard depths
############################################################

fc_spline <- ea_spline(
  obj = fc_spline_input,
  var.name = "THETA",
  d = c(0, 20, 30, 60, 100, 200),
  lam = 0.1,
  vlow = 0,
  show.progress = FALSE
)

############################################################
# Extract 0–20 cm field capacity
############################################################

fc_harmonised <- fc_spline$harmonised

fc_depth020 <- fc_harmonised[, c(
  "id",
  "0-20 cm"
)]

colnames(fc_depth020)[1] <- "PROFILE_ID.x"

############################################################
# Extract profile coordinates
############################################################

fc_coordinates <- fc_profile_coordinates[, c(
  "PROFILE_ID.x",
  "X_WGS84",
  "Y_WGS84"
)]

############################################################
# Merge coordinates with harmonized values
############################################################

fc_profiles <- merge(
  fc_depth020,
  fc_coordinates,
  by = "PROFILE_ID.x"
)

############################################################
# Keep one point per profile
############################################################

fc_profiles_unique <- fc_profiles %>%
  distinct(PROFILE_ID.x, .keep_all = TRUE)
######################################################


############################################################
# Load Field Capacity prediction raster
############################################################

fc_raster <- rast(
  "C:/Users/surya/Downloads/Panos_SHP_data/Resample_1km/Resample_FC_33kpa_PTF07_modified_6_03_tif.tif"
)

############################################################
# Create spatial points
############################################################

fc_points <- vect(
  fc_profiles_unique,
  geom = c("X_WGS84", "Y_WGS84"),
  crs = "EPSG:4326"
)

############################################################
# Extract predicted Field Capacity values
############################################################

fc_prediction <- extract(
  fc_raster,
  fc_points
)

fc_profiles_unique$Predicted_FC <- fc_prediction[, 2]

############################################################
# Remove missing values
############################################################

fc_validation <- fc_profiles_unique %>%
  filter(
    complete.cases(.)
  )

############################################################
# Model performance
############################################################

MAE <- mean(
  abs(
    fc_validation$`0-20 cm` -
      fc_validation$Predicted_FC
  )
)

RMSE <- sqrt(
  mean(
    (
      fc_validation$`0-20 cm` -
        fc_validation$Predicted_FC
    )^2
  )
)

Bias <- mean(
  fc_validation$Predicted_FC -
    fc_validation$`0-20 cm`
)

############################################################
# Print statistics
############################################################
cat("MAE  =", round(MAE, 3), "\n")
cat("RMSE =", round(RMSE, 3), "\n")
cat("Bias =", round(Bias, 3), "\n")

############################################################
# Load EU-SoilHydroGrids raster
############################################################

fc_euhydro_raster <- rast(
  "X:/Users/Surya/EU_SoilHydroGrids_1km/FC_0_20cm_mask.tif"
)

############################################################
# Create spatial points
############################################################

fc_points <- vect(
  fc_validation,
  geom = c("X_WGS84", "Y_WGS84"),
  crs = "EPSG:4326"
)

############################################################
# Extract EU-SoilHydroGrids predictions
############################################################

fc_euhydro_prediction <- extract(
  fc_euhydro_raster,
  fc_points
)

fc_validation$EUHydroGrid_FC <-
  fc_euhydro_prediction[, 2] / 100

############################################################
# Remove missing values
############################################################

fc_validation <- fc_validation %>%
  filter(
    complete.cases(.)
  )

############################################################
# Model performance
############################################################

MAE <- mean(
  abs(
    fc_validation$`0-20 cm` -
      fc_validation$EUHydroGrid_FC
  )
)

RMSE <- sqrt(
  mean(
    (
      fc_validation$`0-20 cm` -
        fc_validation$EUHydroGrid_FC
    )^2
  )
)

Bias <- mean(
  fc_validation$EUHydroGrid_FC -
    fc_validation$`0-20 cm`
)

############################################################
# Print statistics
############################################################

cat("MAE  =", round(MAE, 3), "\n")
cat("RMSE =", round(RMSE, 3), "\n")
cat("Bias =", round(Bias, 3), "\n")


#####################################################

plot_data <- bind_rows(
  data.frame(
    Value = fc_validation$`0-20 cm`,
    Type = "Measured",
    Map = "EUPTFv2–LUCAS map"
  ),
  data.frame(
    Value = fc_validation$Predicted_FC,
    Type = "Predicted",
    Map = "EUPTFv2–LUCAS map"
  ),
  data.frame(
    Value = fc_validation$`0-20 cm`,
    Type = "Measured",
    Map = "EU-SoilHydroGrids map"
  ),
  data.frame(
    Value = fc_validation$EUHydroGrid_FC,
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
    x = expression(FC~"["*cm^3/cm^3*"]"),
    y = "Density"
  )
