##############################################################################
### 03 --- Harmonizing soil data of simulations
##############################################################################
### In this script we harmonize the soil data of the simulations and extract 
### data from Pan-European datasets whenever necessary
### M. Grünig, 06.02.2024
##############################################################################



library(raster)
library(terra)
library(RColorBrewer)
library(sf)
library(dplyr)
library(DBI)
library(stars)
library(ggplot2)
library(exactextractr)



### functions -----------------------------------------------------------------

## WHC calculation
# W. Rammer, April 2019

# Calculate water holding capacity following the approach of Schwalm & Ek (2004, Ecol. Mod.)
# pct_sand, pct_silt, pct_clay: soil texture fractions in %
# soil_depth: stone-free effective soil depth (mm)
# fc_threshold: field capacity (J kg-1 or kPa), default: 15kPa
# pwp_threshold: permanent wilting point (J kg-1 or kPa), default 1500kPa
# return: water holding capacity (between fc and pwp) in mm
calcWHC <- function(pct_sand, pct_silt, pct_clay, soil_depth, fc_threshold=15, pwp_threshold=4000) {
  
  theta_sat = 0.01 * (50.5 - 0.142*pct_sand - 0.037*pct_clay); # Eq. 78
  bt <- 11.43 - 0.103*pct_sand - 0.0687*pct_silt # Eq. 79
  rho_e <- -5.99 + 0.0544*pct_sand + 0.0451*pct_silt # Eq 80
  
  
  fc <- theta_sat * (rho_e / -fc_threshold)^(1/bt) * soil_depth # Eq 76
  pwp <- theta_sat * (rho_e / -pwp_threshold)^(1/bt) * soil_depth # Eq 77
  whc <- fc-pwp
  whc
}



# set the path to directory
public_pth <- "your/path/"

# load europe shapefiles
eu.ext <- c(-10,50,35,80)
eu_poly <- vect(paste0(public_pth, "Projects/Resonate/3pg/gis/europe.shp"))
eu_poly_lr <- vect(paste0(public_pth, "Projects/Resonate/3pg/gis/europe_lowres.shp"))

# projections
proj_forest <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj <- "+proj=longlat +datum=WGS84 +no_defs"


### load data ---

# soil depth
depth <- rast(paste0(public_pth, "/Data/GIS_data/europe/soil/esdac/STU/STU_EU_DEPTH_ROOTS.rst"))
crs(depth) <- proj_forest
depth <- terra::project(depth, proj_forest)
depth <- terra::mask(depth, eu_poly_lr)
depth <- terra::project(depth, proj)
depth[depth == 0] <- NA
depth <- depth * 10
depth_df <- as.data.frame(depth, xy = T)

# clay
clay_s <- rast(paste0(public_pth, "/Data/GIS_data/europe/soil/esdac/STU/STU_EU_T_CLAY.rst"))
crs(clay_s) <- proj_forest
clay_s <- terra::project(clay_s, proj_forest)
clay_s <- terra::mask(clay_s, eu_poly_lr)

# silt
silt_s <- rast(paste0(public_pth, "/Data/GIS_data/europe/soil/esdac/STU/STU_EU_T_SILT.rst"))
crs(silt_s) <- proj_forest
silt_s <- terra::project(silt_s, proj_forest)
silt_s <- terra::mask(silt_s, eu_poly_lr)

# sand
sand_s <- rast(paste0(public_pth, "/Data/GIS_data/europe/soil/esdac/STU/STU_EU_T_SAND.rst"))
crs(sand_s) <- proj_forest
sand_s <- terra::project(sand_s, proj_forest)
sand_s <- terra::mask(sand_s, eu_poly_lr)

# convert to proj
sand_wgs <- terra::project(sand_s, proj)
sand_wgs[sand_wgs == 0] <- NA
silt_wgs <- terra::project(silt_s, proj)
silt_wgs[silt_wgs == 0] <- NA
clay_wgs <- terra::project(clay_s, proj)
clay_wgs[clay_wgs == 0] <- NA
rm(sand_s, silt_s, clay_s)


# water holding capacity from subsoil and topsoil
twhc_s <- rast(paste0(public_pth, "/Data/GIS_data/europe/soil/esdac/STU/STU_EU_S_TAWC.rst"))
twhc_t <- rast(paste0(public_pth, "/Data/GIS_data/europe/soil/esdac/STU/STU_EU_T_TAWC.rst"))
whc <- twhc_s + twhc_t

crs(whc) <- proj_forest
whc <- terra::project(whc, proj_forest)
r_whc <- terra::mask(whc, eu_poly_lr)
r_whc[r_whc == 0] <- NA
r_whc_wgs <- terra::project(r_whc, proj)
rm(r_whc)
whc_df <- as.data.frame(r_whc_wgs, xy = T)


# load texture (if the calculation step was skipped above)
r_texture <- rast(paste0(public_pth, "/Data/GIS_data/europe/soil/texture_categories.tif"))
plot(r_texture)
crs(r_texture) <- proj_forest
r_texture <- terra::project(r_texture, proj_forest, method = "mode")
r_texture[r_texture == 0] <- NA
r_texture <- terra::as.factor(r_texture)


# load nitrogen layer
r_nitro <- rast(paste0(public_pth, "Projects/Resonate/3pg/soil_layers/predicted_N_available_1km_v2_interpolated.tif"))
r_nitro <- terra::project(r_nitro, r_texture)
r_nitro <- terra::mask(r_nitro, r_texture)
r_nitro[r_nitro > 120] <- NA
nitro_wgs <- terra::project(r_nitro, proj)
rm(r_nitro)


### load the simulation metadata -----------------------------------
simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(public_pth, "/Projects/Resonate/data_portal/simulation_data/forest_simulation_db_v1.sqlite"))

tables_con <- dbListTables(simulation_db)
unique(tables_con)

metadata <- dbReadTable(simulation_db, "metadata_all")
dim(metadata)
head(metadata)

# some users used cm instead of mm
cm_users <- c(1003, 1021, 1023, 1032, 1033, 1034, 1035, 1036)
metadata <- metadata %>% 
  mutate(SoilDepth = ifelse(uniqueID %in% cm_users, SoilDepth * 10, SoilDepth))


### Go through data and convert ratings to values
### extract information from pan-European data if needed


# if there is no texture and no soildepth in the data we extract it from the maps
if(nrow(metadata[is.na(metadata$TextureSand),]) > 0){
  
  metadata <- as.data.frame(cbind(metadata, 
                                  sand_extract = terra::extract(sand_wgs, metadata[, c("Lon", "Lat")])[,2],
                                  silt_extract = terra::extract(silt_wgs, metadata[, c("Lon", "Lat")])[,2],
                                  clay_extract = terra::extract(clay_wgs, metadata[, c("Lon", "Lat")])[,2],
                                  depth_extract = terra::extract(depth, metadata[, c("Lon", "Lat")])[,2]))
  
  metadata <- metadata %>% mutate(TextureSand = ifelse(is.na(TextureSand), sand_extract, TextureSand),
                                  TextureSilt = ifelse(is.na(TextureSilt), silt_extract, TextureSilt),
                                  TextureClay = ifelse(is.na(TextureClay), clay_extract, TextureClay),
                                  SoilDepth = ifelse(is.na(SoilDepth), depth_extract, SoilDepth))
  
  
  if(nrow(metadata[is.na(metadata$SoilDepth), ]) > 0){
    
    nas <- metadata[is.na(metadata$SoilDepth), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, depth_buffer = exact_extract(depth, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$SoilDepth), "SoilDepth"] <- as.numeric(nas[, "depth_buffer"])
    
  }
  
  
  if(nrow(metadata[is.na(metadata$SoilDepth), ]) > 0){
    
    nas <- metadata[is.na(metadata$SoilDepth), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 3000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, depth_buffer2 = exact_extract(depth, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$SoilDepth), "SoilDepth"] <- as.numeric(nas[, "depth_buffer2"])
    
  }
  
  
  # same for the other textures
  if(nrow(metadata[is.na(metadata$TextureSand), ]) > 0){
    
    nas <- metadata[is.na(metadata$TextureSand), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, 
                               sand_buffer = exact_extract(sand_wgs, sp_buffer, fun = "mean", progress = F) ,
                               silt_buffer = exact_extract(silt_wgs, sp_buffer, fun = "mean", progress = F) ,
                               clay_buffer = exact_extract(clay_wgs, sp_buffer, fun = "mean", progress = F) ))
    
    metadata[is.na(metadata$TextureSand), "TextureSand"] <- as.numeric(nas[, "sand_buffer"])
    metadata[is.na(metadata$TextureSilt), "TextureSilt"] <- as.numeric(nas[, "silt_buffer"])
    metadata[is.na(metadata$TextureClay), "TextureClay"] <- as.numeric(nas[, "clay_buffer"])
    
  }
  
  if(nrow(metadata[is.na(metadata$TextureSand), ]) > 0){
    
    nas <- metadata[is.na(metadata$TextureSand), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 3000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, 
                               sand_buffer2 = exact_extract(sand_wgs, sp_buffer, fun = "mean", progress = F) ,
                               silt_buffer2 = exact_extract(silt_wgs, sp_buffer, fun = "mean", progress = F) ,
                               clay_buffer2 = exact_extract(clay_wgs, sp_buffer, fun = "mean", progress = F) ))
    
    metadata[is.na(metadata$TextureSand), "TextureSand"] <- as.numeric(nas[, "sand_buffer2"])
    metadata[is.na(metadata$TextureSilt), "TextureSilt"] <- as.numeric(nas[, "silt_buffer2"])
    metadata[is.na(metadata$TextureClay), "TextureClay"] <- as.numeric(nas[, "clay_buffer2"])
    
  }
  
  
}

# check which simulations have no data for WHC but a rating
# rating is from WHC 30mm to 200mm
if(nrow(metadata[is.na(metadata$WHC),]) > 0){
  
  metadata <- metadata %>%
    mutate(WHC = ifelse(is.na(WHC), SoilWaterRating*200, WHC))
  
  metadata <- metadata %>%
    mutate(WHC = ifelse(is.na(WHC),
                        calcWHC(pct_sand = TextureSand,
                                pct_silt = TextureSilt,
                                pct_clay = TextureClay,
                                soil_depth = SoilDepth),
                        WHC))
  
  metadata <- as.data.frame(cbind(metadata, whc_extract = terra::extract(r_whc_wgs, metadata[, c("Lon", "Lat")])[,2]))
  
  metadata <- metadata %>%
    mutate(WHC = ifelse(is.na(WHC) & is.na(SoilDepth), whc_extract, WHC))
  
  # if there are mssing points - get from raster layers
  # there is a handful (12) simulations in two locations that are on the coast in northern germany
  # and get a NA from the extract function, we make correction by adding a 1km buffer around the locations
  
  
  if(nrow(metadata[is.na(metadata$WHC), ]) > 0){
    
    nas <- metadata[is.na(metadata$WHC), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, whc_buffer = exact_extract(r_whc_wgs, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$WHC), "WHC"] <- as.numeric(nas[, "whc_buffer"])
    
  }
}



# if available nitrogen is missing we take the fertiliy rating
# for some cases there is no rating, therefore we extract values from the map
if(nrow(metadata[is.na(metadata$AvailableNitrogen), ]) > 0){
  
  # transform the rating value to what we assigned them in the model data interface
  metadata <- metadata %>%
    mutate(AvailableNitrogen = ifelse(is.na(AvailableNitrogen), 
                                      ifelse(FertilityRating < 0.5, FertilityRating*80 + 20*(1-FertilityRating),
                                             FertilityRating*100), 
                                      AvailableNitrogen))
  
  metadata <- as.data.frame(cbind(metadata, 
                                  nitro_extract = terra::extract(nitro_wgs, metadata[, c("Lon", "Lat")])[,2]))
  
  metadata <- metadata %>% 
    mutate(AvailableNitrogen = ifelse(is.na(AvailableNitrogen), nitro_extract, AvailableNitrogen))
  
  # there is a handful (12) simulations in two locations that are on the coast in northern germany
  # and get a NA from the extract function, we make correction by adding a 1km buffer around the locations
  
  if(nrow(metadata[is.na(metadata$AvailableNitrogen), ]) > 0){
    
    nas <- metadata[is.na(metadata$AvailableNitrogen), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, nitro_buffer = exact_extract(nitro_wgs, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$AvailableNitrogen), "AvailableNitrogen"] <- as.numeric(nas[, "nitro_buffer"])
    
  }
  
  if(nrow(metadata[is.na(metadata$AvailableNitrogen), ]) > 0){
    
    nas <- metadata[is.na(metadata$AvailableNitrogen), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1500) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, nitro_buffer2 = exact_extract(nitro_wgs, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$AvailableNitrogen), "AvailableNitrogen"] <- as.numeric(nas[, "nitro_buffer2"])
    
  }
  
}





head(metadata)

# check that no more NAs are in the dataset for the needed cols
nrow(metadata[is.na(metadata$AvailableNitrogen) | is.na(metadata$WHC) | is.na(metadata$SoilDepth) | is.na(metadata$TextureSand), ])

metadata <- metadata[, 1:19]
dbWriteTable(conn = simulation_db, name = "metadata_harmonized", value = metadata, overwrite = T)
dbDisconnect(simulation_db)

