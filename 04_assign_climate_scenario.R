##################################################
### 04 --- get best fitting climate scenario
##################################################
### In this script we obtain the best fitting climate 
### scenario from the climate DB for each simulation
### 
### M. Grünig, 06.02.2024
##################################################


library(raster)
library(terra)
library(RColorBrewer)
library(rgdal)
library(sf)
library(dplyr)
library(DBI)
library(stars)
library(ggplot2)
library(exactextractr)
library(sparklyr)
library(arrow)
library(data.table)
library(collapse)
library(parallel)


### load reference grid
source("functions/get_point_id.R")

path <- "/your/path/"


### load the simulation DB -----------------------------------
simulation_db <- DBI::dbConnect(RSQLite::SQLite(),
                                paste0(path, "/forest_simulation_db_v1.sqlite"))
tables_con <- dbListTables(simulation_db)

# load metadata and check
metadata <- dbReadTable(simulation_db, "metadata_harmonized")
dim(metadata)
head(metadata)


# load the climate slopes, extracted from the climate database in the previous script
climate_slopes <- dbReadTable(simulation_db, "climate_slopes_coordinates") %>% data.table()

# isolate the simulation data tables
sim_tables <- tables_con[grepl(paste0("simulation"), tables_con)]

# create rcp and gcm list to choose from
rcp_list <- list(c("historical"), c("rcp_2_6"), c("rcp_4_5"), c("rcp_8_5"),
                 c("rcp_4_5", "rcp_8_5"), c("rcp_2_6", "rcp_4_5", "rcp_8_5"))

gcm_list <- list(c("ICHEC-EC-EARTH"), c("MPI-M-MPI-ESM-LR"), c("NCC-NorESM1-M"),
                 c("ICHEC-EC-EARTH", "MPI-M-MPI-ESM-LR", "NCC-NorESM1-M")) 



### Create a function that compares the slope of the simulation with slopes from the different climate scenarios -------

scenario_compare <- function(x, column_names){
  
  x <- as.vector(x)
  names(x) <- column_names
  sim_id <- as.numeric(x["ID"])
  
  # get the climate data from the simulation
  simdata <- simulation_data %>%
    filter(simulationID == as.numeric(sim_id)) %>% collect()
  
  # if empty return NAs
  if(nrow(simdata) == 0){return(c(NA, NA, NA))}
  
  # get the temp and precip and calibrate linear models
  simdata <- simdata %>% fsubset(!is.na(Temp) & !is.na(Precip))
  sim_temp_mod <- flm(Temp ~ Year, data = simdata)
  sim_prec_mod <- flm(Precip ~ Year, data = simdata)
  
  # calculate averages
  sim_temp_avg <- mean(simdata$Temp)
  sim_prec_avg <- mean(simdata$Precip)
  
  # get climate for the gridcell
  sim_coords <- x[c("Lon", "Lat")]
  grid_id <- get_point_id(sim_coords)
  
  # check which rcp we have to look at
  rcp_indicator <- ifelse(grepl("aseline", x[c("Climate")]) | grepl("bserved", x[c("Climate")]), 1, 
                          ifelse(grepl("26", x[c("Climate")]) | grepl("2.6", x[c("Climate")]), 2, 
                                 ifelse(grepl("45", x[c("Climate")]) | grepl("4.5", x[c("Climate")]), 3, 
                                        ifelse(grepl("85", x[c("Climate")]) | grepl("8.5", x[c("Climate")]), 4,
                                               ifelse(grepl("6.0", x[c("Climate")]) | grepl("60", x[c("Climate")]), 5,
                                                      ifelse(grepl("SRES", x[c("Climate")]), 5, 6))))))
  
  rcps <- rcp_list[[rcp_indicator]]
  
  # check if there is information on the gcm in the metadata
  gcm_indicator <- ifelse(is.na(x[c("GCM")]), 4,
                          ifelse(grepl("ICHEC", x[c("GCM")]) | grepl("EARTH", x[c("GCM")]) | grepl("earth", x[c("GCM")]), 1, 
                                 ifelse(grepl("MPI", x[c("GCM")]) | grepl("mpi", x[c("GCM")]), 2, 
                                        ifelse(grepl("NCC", x[c("GCM")]) | grepl("NorESM1", x[c("GCM")]) | grepl("ncc", x[c("GCM")]), 3, 4))))
  
  gcms <- gcm_list[[gcm_indicator]]
  
  # create empty output DF
  diff_df <- NULL
  
  # loop over the potential gcms and rcps
  for(g in gcms){
    
    for(c in rcps){
      
      # for the coordinates of the simulation we get the precip and temp slope from the climate database for the focal scenario
      focal_gridcell <- climate_slopes %>% collapse::fsubset(Lon == as.numeric(sim_coords[1])) %>% collapse::fsubset(Lat == as.numeric(sim_coords[2])) %>% distinct() %>% as.data.frame()
      temp_slope_clim <- as.numeric(focal_gridcell[colnames(focal_gridcell) == paste0("temp_slope_", gsub("-", ".", g), "_", c)])
      prec_slope_clim <- as.numeric(focal_gridcell[colnames(focal_gridcell) == paste0("prec_slope_", gsub("-", ".", g), "_", c)])
      
      # calculate the deltas. temp is additive, precip multiplicative
      delta_t <- abs(as.numeric(sim_temp_mod[2,1]) - temp_slope_clim)
      #delta_p <- abs(1 - as.numeric(sim_prec_mod[2,1]) / prec_slope_clim) # before we had as.numeric(sim_prec_mod[2,1]) / prec_slope_clim but this is not correct I think
      # the alternative is to do it like for 
      delta_p <- abs(as.numeric(sim_prec_mod[2,1]) - prec_slope_clim)
      
      # add together with a arbitrary 0.1 weight for precip
      dev_metric <- delta_t + 0.1 * delta_p
      dev_metric <- data.table(dev_metric)
      colnames(dev_metric) <- paste0(g, "_", c)
      
      diff_df <- cbind(diff_df, dev_metric)
    }
  }
  
  # get the closest scenario
  min.val <- colnames(diff_df)[which.min(diff_df)]
  
  # calculate the offset from the averages
  temp_offset <- as.numeric(focal_gridcell[colnames(focal_gridcell) == paste0("temp_avg_", gsub("-", ".", min.val))]) - sim_temp_avg
  prec_offset <- as.numeric(focal_gridcell[colnames(focal_gridcell) == paste0("prec_avg_", gsub("-", ".", min.val))])/sim_prec_avg
  
  # return the values for scenario, temp offset and prec offset
  if(is.na(min.val)){return(c(NA, temp_offset, prec_offset))}else{
    return(c(as.character(min.val), as.numeric(temp_offset), as.numeric(prec_offset)))
  }
  
  gc()
}

# compile the function to make it faster
library(compiler)
scenario_compare_cmp <- cmpfun(scenario_compare)


### run the function ----------------------------------------------------------------------------------
usernames <- read.csv(paste0(path, "/user_ids.csv"), sep=";")

# list all users
users <- list.files(paste0(path))
target_scen_all_list <- list()

system.time({
  
  r <- 0  
  for(u in users){
    if(u %in% empty_users){next}
    
    print(u)
    r <- r + 1
    
    #open database
    simulation_db <- DBI::dbConnect(RSQLite::SQLite(),
                                    paste0(path, "/forest_simulation_db_v1.sqlite"))
    
    userID <- gsub("user_", "", u)
    
    # subset meta data
    metadata_sub <- metadata %>% collapse::fsubset(uniqueID == userID)
    
    # if there is nothing, we skip
    if(nrow(metadata_sub) == 0){
      cat(paste0(u))
      next}
    
    # convert to list to enable mclapply
    metadata_sub_list <- split(metadata_sub, seq(nrow(metadata_sub)))  
    
    # load the simulation data of the user
    simulation_data <- tbl(simulation_db, sim_tables[grepl(userID, sim_tables)]) %>%
      dplyr::select(simulationID, Year, Temp, Precip) %>% collect()
    
    # apply function in parallel
    target_scen <- parallel::mclapply(metadata_sub_list, FUN = scenario_compare_cmp,
                                      column_names = names(metadata_sub_list[[1]]), mc.cores = 8)
    
    # combine to a table 
    target_scen_df <- data.table(do.call(rbind, target_scen))
    colnames(target_scen_df) <- c("scenario", "temp_offset", "prec_offset")
    target_scen_df <- target_scen_df %>% mutate_at(c('temp_offset', 'prec_offset'), as.numeric)
    
    # add to list of users  
    target_scen_all_list[[r]] <- target_scen_df
    
    # close DB  
    dbDisconnect(simulation_db)
    
  }# close users loop
  
})

# combine all tables
target_scen_all <- data.table(do.call(rbind, target_scen_all_list))
dim(target_scen_all)

# bind together to metadata and save in the SQLite DB
metadata_fit_scenario <- cbind(metadata, target_scen_all)
metadata_fit_scenario <- metadata_fit_scenario %>% filter(!is.na(scenario))
head(metadata_fit_scenario)

simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.sqlite"))
dbWriteTable(conn = simulation_db, name = "metadata_climate_scenario",
             value = metadata_fit_scenario, overwrite = T)

dbDisconnect(simulation_db)



