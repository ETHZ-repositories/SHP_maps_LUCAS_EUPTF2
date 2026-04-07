################PTF_making_maps_Europe_ Uisng soil texture, bulk density, depth and soil carbon maps
Thetas_PTF<- THS_PTF03

Thetas_PTF$variable.importance

grid <- list.files("C:/Users/surya/Downloads/SHP_data/" , pattern = "*.tif$")

All_cov <- raster::stack(paste0("C:/Users/surya/Downloads/SHP_data/", grid))

names(All_cov) <- tools::file_path_sans_ext(grid)

p2 = predict(All_cov,Thetas_PTF, progress='window',type = "response",fun = function(model, ...) predict(model, ...)$predictions)

writeRaster(p2, "C:/Users/surya/Downloads/Panos_SHP_data/Thetas_PTF03.tif")

p7_lower = predict(All_cov,Thetas_PTF, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.05)$predictions)

p7_upper = predict(All_cov,Thetas_PTF, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.95)$predictions)

plot(p7_upper)

plot(p7_lower)

Uncertainity_Thetas_PTF<- p7_upper - p7_lower

plot(Uncertainity_Thetas_PTF)

writeRaster(Uncertainity_Thetas_PTF, "C:/Users/surya/Downloads/Panos_SHP_data/Folder_data/Maps_withOC/Uncertainity_Thetas_PTF_modified_6_03.tif")

##################################################### with OC PTFs

FC_PTF07_33kpa<- FC_PTF07

FC_PTF07_33kpa$variable.importance

grid <- list.files("C:/Users/surya/Downloads/SHP_data/Modified_files/" , pattern = "*.tif$")

All_cov <- raster::stack(paste0("C:/Users/surya/Downloads/SHP_data/Modified_files/", grid))

names(All_cov) <- tools::file_path_sans_ext(grid)

p4_33_07 = predict(All_cov,FC_PTF07_33kpa, progress='window',type = "response",fun = function(model, ...) predict(model, ...)$predictions)

writeRaster(p4_33_07, "C:/Users/surya/Downloads/SHP_data/Folder_data/Maps_withOC/FC_33kpa_PTF07_modified_6_03.tif")

p3_lower = predict(All_cov,FC_PTF07_33kpa, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.05)$predictions)

p3_upper = predict(All_cov,FC_PTF07_33kpa, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.95)$predictions)

plot(p3_upper)

plot(p3_lower)

Uncertainity_FC_PTF07_33kpa<- p3_upper - p3_lower

plot(Uncertainity_FC_PTF07_33kpa)

writeRaster(Uncertainity_FC_PTF07_33kpa, "C:/Users/surya/Downloads/SHP_data/Folder_data/Maps_withOC/Uncertainity_FC_33kpa_PTF07_modified_6_03_C.tif")

#######################################################

FC_PTF07_10kpa<- FC2_PTF07

FC_PTF07_10kpa$variable.importance

grid <- list.files("C:/Users/surya/Downloads/SHP_data/Modified_files/" , pattern = "*.tif$")

All_cov <- raster::stack(paste0("C:/Users/surya/Downloads/SHP_data/Modified_files/", grid))

names(All_cov) <- tools::file_path_sans_ext(grid)

p4_10_07 = predict(All_cov,FC_PTF07_10kpa, progress='window',type = "response",fun = function(model, ...) predict(model, ...)$predictions)

writeRaster(p4_10_07, "C:/Users/surya/Downloads/SHP_data/Folder_data/Maps_withOC/FC_10kpa_PTF07_modified_6_03.tif")

p4_lower = predict(All_cov,FC_PTF07_10kpa, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.05)$predictions)

p4_upper = predict(All_cov,FC_PTF07_10kpa, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.95)$predictions)

plot(p4_upper)

plot(p4_lower)

Uncertainity_FC_PTF07_10kpa<- p4_upper - p4_lower

plot(Uncertainity_FC_PTF07_10kpa)

writeRaster(Uncertainity_FC_PTF07_10kpa, "C:/Users/surya/Downloads/SHP_data/Folder_data/Maps_withOC/Uncertainity_FC_10kpa_PTF07_modified_6_03.tif")



#######################################################

WP_PTF07_150kpa<- WP_PTF07

WP_PTF07_150kpa$variable.importance

grid <- list.files("C:/Users/surya/Downloads/SHP_data/Modified_files/" , pattern = "*.tif$")

All_cov <- raster::stack(paste0("C:/Users/surya/Downloads/SHP_data/Modified_files/", grid))

names(All_cov) <- tools::file_path_sans_ext(grid)

p4_150_07 = predict(All_cov,WP_PTF07_150kpa, progress='window',type = "response",fun = function(model, ...) predict(model, ...)$predictions)

writeRaster(p4_150_07, "C:/Users/surya/Downloads/SHP_data/Folder_data/Maps_withOC/WP_150kpa_PTF07_modified_6_03.tif")

plot(p4_150_07)

p5_lower = predict(All_cov,WP_PTF07_150kpa, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.05)$predictions)

p5_upper = predict(All_cov,WP_PTF07_150kpa, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.95)$predictions)

plot(p5_upper)

plot(p5_lower)

Uncertainity_WP_PTF07_150kpa<- p5_upper - p5_lower

plot(Uncertainity_WP_PTF07_150kpa)

writeRaster(Uncertainity_WP_PTF07_150kpa, "C:/Users/surya/Downloads/SHP_data/Folder_data/Maps_withOC/Uncertainity_WP_PTF07_150kpa_modified_6_03.tif")


#######################################################

Ks_PTF02<- KS_PTF02

Ks_PTF02$variable.importance

grid <- list.files("C:/Users/surya/Downloads/SHP_data/Layer_based_OC/Ks_layer/" , pattern = "*.tif$")

All_cov <- raster::stack(paste0("C:/Users/surya/Downloads/SHP_data/Layer_based_OC/Ks_layer/", grid))

names(All_cov) <- tools::file_path_sans_ext(grid)

p4_Ks_02 = predict(All_cov,Ks_PTF02, progress='window',type = "response",fun = function(model, ...) predict(model, ...)$predictions)

writeRaster(p4_Ks_02, "C:/Users/surya/Downloads/SHP_data/Folder_data/Maps_withOC/Ks_PTF02_modified_6_03.tif")

plot(p4_Ks_02)

p6_lower = predict(All_cov,Ks_PTF02, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.05)$predictions)

p6_upper = predict(All_cov,Ks_PTF02, progress='window',type = "quantiles",fun = function(model, ...) predict(model, ..., quantiles = 0.95)$predictions)

plot(p6_upper)

plot(p6_lower)

Uncertainity_Ks_PTF02<- p6_upper - p6_lower

plot(Uncertainity_Ks_PTF02)

writeRaster(Uncertainity_Ks_PTF02, "C:/Users/surya/Downloads/SHP_data/Folder_data/Maps_withOC/Uncertainity_Ks_PTF02_modified_6_03.tif")


# ################################################PTF01
# 
# FC_PTF01_33<- FC_PTF01
# 
# FC_PTF01_33$variable.importance
# 
# grid <- list.files("C:/Users/surya/Downloads/Panos_SHP_data/PSD_covariates/" , pattern = "*.tif$")
# 
# All_cov <- raster::stack(paste0("C:/Users/surya/Downloads/Panos_SHP_data/PSD_covariates/", grid))
# 
# names(All_cov) <- tools::file_path_sans_ext(grid)
# 
# p4_PTF_01 = predict(All_cov,FC_PTF01_33, progress='window',type = "response",fun = function(model, ...) predict(model, ...)$predictions)
# 
# writeRaster(p4_PTF_01, "C:/Users/surya/Downloads/Panos_SHP_data/PSD_covariates/FC_PTF01_33kpa.tif")
# 
# plot(p4_PTF_01)
