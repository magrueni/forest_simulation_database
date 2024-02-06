##################################################
### 05 --- extract the timeseries for simulations
##################################################
### In this script we extract the matching time-
### series for each simulation run

### M. Grünig, 06.02.2024
##################################################


library(raster)
library(terra)
library(sf)
library(dplyr)
library(DBI)
library(stars)
library(ggplot2)
library(exactextractr)
library(arrow)
library(data.table)
library(collapse)
library(future.apply)
library(compiler)
library(parallel)
library(dbplyr)


path <- "/your/path/"



### load functions ---------------------------------------------------------------------

# get pointid func
source("functions/get_point_id.R")


# create time series function. This function compares the temp and prec from the 
# scenarios with the one from the simulations
time_series_creation <- function(x, temp, prec){
  
  dev_metric <- unlist(as.vector(abs(as.numeric(x["Temp"]) - temp) + (0.1 * as.numeric(x["Precip"])/prec)))
  
  min.year.val <- sample(names(base::sort(dev_metric)[1:3]), 1)
  min.year.val <- as.numeric(gsub("X", "", min.year.val, perl = TRUE))
  
  return(min.year.val)
  
}

time_series_creation_cmp <- cmpfun(time_series_creation)


# core function
# this function takes the climate data of each simulation and 
# uses then the function above to match it with the best fitting scenario
core_function <- function(meta_data, col_names){
  
  # obtain the sim id
  meta_data <- as.vector(meta_data)
  names(meta_data) <- col_names
  sim_id <- as.numeric(meta_data["ID"])
  
  # get the simulation data and extract temp and precip
  simdata_sub <- simulation_data %>%
    filter(simulationID == as.numeric(sim_id)) %>%
    collapse::fsubset(!is.na(Temp) & !is.na(Precip)) %>%
    collapse::fselect(Year, Temp, Precip)
  
  # get scenrio grid id
  sim_coords <- as.numeric(meta_data[c("Lon", "Lat")])
  grid_id <- get_point_id(sim_coords)
  
  # get the information of the scenario from the metadata
  target_scen <- as.character(meta_data["scenario"])
  target_scen <- gsub("-", "_", target_scen)
  
  temp_data_scen <- eval(parse(text = paste0("temp_", target_scen)))
  prec_data_scen <- eval(parse(text = paste0("prec_", target_scen)))
  
  # subset the temp data for the coordinate
  temp_data <- temp_data_scen %>% 
    collapse::fsubset(point_id == grid_id) 
  temp_data <- data.table(unique(temp_data[, 4:ncol(temp_data)]))
  temp_data <- temp_data - meta_data$temp_offset
  
  # subset the prec data 
  prec_data <- prec_data_scen %>% 
    collapse::fsubset(point_id == grid_id) 
  prec_data <- data.table(unique(prec_data[, 4:ncol(prec_data)]))
  prec_data <- prec_data / meta_data$prec_offset
  
  # apply matching function
  new_time_series <- apply(simdata_sub, 1, time_series_creation_cmp, temp = temp_data, prec = prec_data)
  
  return(new_time_series)

}

myFuncCmp <- cmpfun(core_function)


### load data from sqlite -----------------------------------------------------------------------------------

# open connection
simulation_db <- DBI::dbConnect(RSQLite::SQLite(),
                                paste0(path, "/forest_simulation_db_v1.sqlite"))

# read tables
tables_con <- dbListTables(simulation_db)

# load metadata
metadata <- dbReadTable(simulation_db, "metadata_climate_scenario")

# load climate slopes that were calculated in script "06_climateDB_slopes.R"
climate_slopes <- dbReadTable(simulation_db, "climate_slopes_coordinates") %>%
  data.table()

# get all simulation data tables
sim_tables <- tables_con[grepl(paste0("simulation"), tables_con)]

# define all scenarios
gcms <- c("ICHEC-EC-EARTH", "MPI-M-MPI-ESM-LR", "NCC-NorESM1-M")
rcps <- c("rcp_2_6", "rcp_4_5", "rcp_8_5", "historical")


# loop over scenarios to read in climate means
for(g in gcms){
  
  for( c in rcps){
    
    # combine gcms and rcps to scneario names
    scen_name <- paste0(g, "_", c)
    scen_name <- gsub("-", "_", scen_name)
    
    # read in climate means
    temp_data <- dbReadTable(simulation_db, tables_con[grepl(paste0("annual_mean_temp_", g, "_", c), tables_con)])
    assign(paste0("temp_", scen_name), temp_data)
    prec_data <- dbReadTable(simulation_db, tables_con[grepl(paste0("annual_precip_", g, "_", c), tables_con)])
    assign(paste0("prec_", scen_name), prec_data)
    
  }
}


# get list of users we need to loop over
usernames <- read.csv(paste0(path, "/user_ids.csv"), sep=";")
path_responses <- paste0(path, "/responses/")

# adapt depending on iteration
users <- list.files(paste0(path_responses))

  
# loop over users
for(u in users){
  if(u %in% empty_users){next}
  
  # connect to DB
  simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.sqlite"))
  
  # monitor progress
  print(u)
  userID <- gsub("user_", "", u)
  
  # subset meta data
  metadata_sub <- metadata %>% collapse::fsubset(uniqueID == userID) %>% filter(!is.na(scenario))
  
  # if empty metadata skip
  if(nrow(metadata_sub) == 0){
    cat(paste0(u))
    next}
  
  # for the big user 
  if(nrow(metadata_sub) > 50000){
    
    chnks <- split(c(1:nrow(metadata_sub)), ceiling(seq_along(c(1:nrow(metadata_sub)))/50000))
    
    # load simulation data
    simulation_data <- tbl(simulation_db, sim_tables[grepl(userID, sim_tables)]) %>%
      dplyr::select(simulationID, Year, Temp, Precip) %>% collect()
    
    
    for(chn in 1:length(chnks)){
      
      print(paste0(chn, " / ", length(chnks)))
      
      metadata_sub_chnk <- metadata_sub[chnks[[chn]][1]:chnks[[chn]][length(chnks[[chn]])], ]
      
      # convert to list
      metadata_sub_list <- split(metadata_sub_chnk, seq(nrow(metadata_sub_chnk)))  
      
      # lapply over all metadata entries and apply matching functions defined above
      set.seed(10)
      new_series <- parallel::mclapply(metadata_sub_list, FUN = myFuncCmp, col_names = names(metadata_sub_list[[1]]), mc.cores = 16)
      
      # rearrange 
      new_series_df <- data.table(t(sapply(new_series, "length<-", max(lengths(new_series)))))
      colnames(new_series_df) <- paste0("Year", seq_along(new_series_df))
      
      # add user ID
      time_series_df <- cbind(uniqueID = rep(userID, nrow(metadata_sub_chnk)), simulationID = metadata_sub_chnk$ID, new_series_df)
      
      # write out to database
      if(chn == 1){
        dbWriteTable(conn = simulation_db, name = paste0(u, "_new_timeseries"), value = time_series_df, overwrite = T)
      }else{
        dbWriteTable(conn = simulation_db, name = paste0(u, "_new_timeseries"), value = time_series_df, overwrite = F, append = T)
      }
      
    }
    
    # clean up
    rm(simulation_data)
    gc()
    
    # disconnect 
    dbDisconnect(simulation_db)
    
  }else{
    
    # convert to list
    metadata_sub_list <- split(metadata_sub, seq(nrow(metadata_sub)))  
    
    # load simulation data
    simulation_data <- tbl(simulation_db, sim_tables[grepl(userID, sim_tables)]) %>%
      dplyr::select(simulationID, Year, Temp, Precip) %>% collect()
    
    # lapply over all metadata entries and apply matching functions defined above
    set.seed(10)
    new_series <- parallel::mclapply(metadata_sub_list, FUN = myFuncCmp, col_names = names(metadata_sub_list[[1]]), mc.cores = 16)
    
    # rearrange 
    new_series_df <- data.table(t(sapply(new_series, "length<-", max(lengths(new_series)))))
    colnames(new_series_df) <- paste0("Year", seq_along(new_series_df))
    
    # add user ID
    time_series_df <- cbind(uniqueID = rep(userID, nrow(metadata_sub)), simulationID = metadata_sub$ID, new_series_df)
    
    # write out to database
    dbWriteTable(conn = simulation_db, name = paste0(u, "_new_timeseries"), value = time_series_df, overwrite = T)
    
    # clean up
    rm(simulation_data)
    gc()
    
    # disconnect 
    dbDisconnect(simulation_db)
    
  }
}
