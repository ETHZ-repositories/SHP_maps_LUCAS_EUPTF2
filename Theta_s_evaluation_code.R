############################################################
# Validation of Saturated Water Content (ThetaS)

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
# Select saturated water content (HEAD = 0 cm)
############################################################

thetas_measurements <- retention_data %>%
  filter(HEAD == 0)

thetas_measurements <- thetas_measurements[, 1:4]

############################################################
# Average duplicate measurements
############################################################

thetas_measurements_mean <- aggregate(
  . ~ SAMPLE_ID,
  data = thetas_measurements,
  mean
)

############################################################
# Merge profile information
############################################################

thetas_profile_data <- merge(
  general_data,
  thetas_measurements_mean,
  by = "PROFILE_ID"
)

thetas_profile_basic <- merge(
  thetas_profile_data,
  basic_data,
  by = "SAMPLE_ID"
)

############################################################
# Keep required columns
############################################################

thetas_profile_coordinates <- thetas_profile_basic[, c(
  1, 2,
  6, 7,
  9, 10,
  13, 14
)]

thetas_profile_coordinates <- thetas_profile_coordinates %>%
  filter(!is.na(X_WGS84))

############################################################
# Prepare spline input
############################################################

thetas_spline_input <- thetas_profile_coordinates[, c(
  "PROFILE_ID.x",
  "SAMPLE_DEP_TOP",
  "SAMPLE_DEP_BOT",
  "THETA"
)]

colnames(thetas_spline_input) <- c(
  "PROFILE_ID",
  "TOP",
  "BOT",
  "THETA"
)

############################################################
# Remove invalid horizons
############################################################

thetas_spline_input <- thetas_spline_input %>%
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

thetas_spline_input <- thetas_spline_input %>%
  group_by(PROFILE_ID, TOP, BOT) %>%
  summarise(
    THETA = mean(THETA, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
# Keep profiles with at least two horizons
############################################################

thetas_spline_input <- thetas_spline_input %>%
  group_by(PROFILE_ID) %>%
  filter(n() >= 2) %>%
  ungroup()

############################################################
# Harmonize to standard depths
############################################################

thetas_spline <- ea_spline(
  obj = thetas_spline_input,
  var.name = "THETA",
  d = c(0, 20, 30, 60, 100, 200),
  lam = 0.1,
  vlow = 0,
  show.progress = FALSE
)

############################################################
# Extract 0–20 cm ThetaS
############################################################

thetas_harmonised <- thetas_spline$harmonised

thetas_depth020 <- thetas_harmonised[, c(
  "id",
  "0-20 cm"
)]

colnames(thetas_depth020)[1] <- "PROFILE_ID.x"

############################################################
# Extract profile coordinates
############################################################

thetas_coordinates <- thetas_profile_coordinates[, c(
  "PROFILE_ID.x",
  "X_WGS84",
  "Y_WGS84"
)]

############################################################
# Merge coordinates with harmonized values
############################################################

thetas_profiles <- merge(
  thetas_depth020,
  thetas_coordinates,
  by = "PROFILE_ID.x"
)

############################################################
# Keep one point per profile
############################################################

thetas_profiles_unique <- thetas_profiles %>%
  distinct(PROFILE_ID.x, .keep_all = TRUE)
######################################################


############################################################
# Load ThetaS prediction raster
############################################################

thetas_raster <- rast(
  "C:/Users/surya/Downloads/Panos_SHP_data/Resample_1km/Resample_Thetas_PTF03_tif.tif"
)

############################################################
# Create spatial points
############################################################

thetas_points <- vect(
  thetas_profiles_unique,
  geom = c("X_WGS84", "Y_WGS84"),
  crs = "EPSG:4326"
)

############################################################
# Extract predicted ThetaS values
############################################################

thetas_prediction <- extract(
  thetas_raster,
  thetas_points
)

thetas_profiles_unique$Predicted_ThetaS <-
  thetas_prediction[, 2]

############################################################
# Remove missing values
############################################################

thetas_validation <- thetas_profiles_unique %>%
  filter(
    complete.cases(.)
  )

############################################################
# Model performance
############################################################

MAE <- mean(
  abs(
    thetas_validation$`0-20 cm` -
      thetas_validation$Predicted_ThetaS
  )
)

RMSE <- sqrt(
  mean(
    (
      thetas_validation$`0-20 cm` -
        thetas_validation$Predicted_ThetaS
    )^2
  )
)

Bias <- mean(
  thetas_validation$Predicted_ThetaS -
    thetas_validation$`0-20 cm`
)

############################################################
# Print statistics
############################################################

cat("MAE  =", round(MAE, 3), "\n")
cat("RMSE =", round(RMSE, 3), "\n")
cat("Bias =", round(Bias, 3), "\n")

###################################

######################################

############################################################
# Load EU-SoilHydroGrids raster
############################################################

thetas_euhydro_raster <- rast(
  "X:/Users/Surya/EU_SoilHydroGrids_1km/Ths_0_20cm_resample_mask.tif"
)

############################################################
# Create spatial points
############################################################

thetas_points <- vect(
  thetas_validation,
  geom = c("X_WGS84", "Y_WGS84"),
  crs = "EPSG:4326"
)

############################################################
# Extract EU-SoilHydroGrids predictions
############################################################

thetas_euhydro_prediction <- extract(
  thetas_euhydro_raster,
  thetas_points
)

thetas_validation$EUHydroGrid_ThetaS <-
  thetas_euhydro_prediction[, 2] / 100

############################################################
# Model performance
############################################################


MAE <- mean(
  abs(
    thetas_validation$`0-20 cm` -
      thetas_validation$EUHydroGrid_ThetaS
  )
)

RMSE <- sqrt(
  mean(
    (
      thetas_validation$`0-20 cm` -
        thetas_validation$EUHydroGrid_ThetaS
    )^2
  )
)

Bias <- mean(
  thetas_validation$EUHydroGrid_ThetaS -
    thetas_validation$`0-20 cm`
)

############################################################
# Print statistics
############################################################


cat("MAE  =", round(MAE, 3), "\n")
cat("RMSE =", round(RMSE, 3), "\n")
cat("Bias =", round(Bias, 3), "\n")

#######################################################


plot_data <- bind_rows(
  data.frame(
    Value = thetas_validation$`0-20 cm`,
    Type = "Measured",
    Map = "EUPTFv2–LUCAS map"
  ),
  data.frame(
    Value = thetas_validation$Predicted_ThetaS,
    Type = "Predicted",
    Map = "EUPTFv2–LUCAS map"
  ),
  data.frame(
    Value = thetas_validation$`0-20 cm`,
    Type = "Measured",
    Map = "EU-SoilHydroGrids map"
  ),
  data.frame(
    Value = thetas_validation$EUHydroGrid_ThetaS,
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
    x = expression(THS~"["*cm^3/cm^3*"]"),
    y = "Density"
  )
