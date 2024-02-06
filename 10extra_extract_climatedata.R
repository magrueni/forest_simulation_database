##################################################
### 08 --- extract the timeseries for simulations
##################################################
### script for extracting climate data for each 
### simulation and compress it

### M. Grünig, 06.02.2024
##################################################


# note that this script does not work without access to the climate database
# which we can't provide. 


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
library(compiler)
library(parallel)

### spark config ---------------------------------------------------------------
# Set memory allocation for whole local Spark instance
Sys.setenv("SPARK_MEM" = "26g")

config <- spark_config()
config$`spark.driver.memory` <- "24g"
config$`sparklyr.shell.driver-memory` <- "24g"
config$`spark.executor.memory` <- "24g"
config$`spark.driver.maxResultSize` <- "12g"
config$`spark.yarn.executor.memoryOverhead` <- "4g"
config$sparklyr.gateway.port = 8890
config$sparklyr.gateway.start.timeout = 180
config$sparklyr.gateway.connect.timeout = 180
config$`sparklyr.cores.local` <- 8
config$`spark.executor.cores` <- 8
config$spark.executor.extraJavaOptions <- "-XX:MaxHeapFreeRatio=60 -XX:MinHeapFreeRatio=30"    # Set garbage collection settings
config$spark.driver.extraJavaOptions <- "-XX:MaxHeapFreeRatio=60 -XX:MinHeapFreeRatio=30"    # Set garbage collection settings

options(java.parameters = "-Xmx8048m")

### functions --------------------
path <- "/your/path/"
x_rast <- rast("/reference_grid.tif")
crs(x_rast) <- "+proj=longlat +datum=WGS84 +no_defs"
ref_tab <- read.csv("/reference_grid_tab.csv")


get_point_id <- function(x){
  
  # extract the grid ID
  lon <- as.numeric(x[1])
  lat <- as.numeric(x[2])
  
  coords <- as.data.frame(cbind(lon, lat))
  
  point_extr <- terra::extract(x_rast, coords)[2]
  point_extr <- as.numeric(point_extr)
  
  if(is.na(point_extr)){
    closest_grid <- ref_tab[which.min(abs(lon - ref_tab$wgs_x) + abs(lat - ref_tab$wgs_y)), ]
    point_extr <- as.numeric(closest_grid[1])}
  
  return(point_extr)
  
}


# climate extraction function
extract_clim_func <- function(foc_year, clim_tab, sim_point, temp_offset, prec_offset){
  
  # clim_tab_sim <- clim_tab %>% fsubset(point_id == sim_point) %>% fsubset(year == foc_year) %>% fselect(tas, prec, rad, vpd)
  # 
  clim_tab_sim <- clim_tab[point_id == sim_point & year == foc_year, .(tas, prec, rad, vpd)]
  ### ADD THE OFFSET!
  clim_tab_sim[, c("tas", "prec") := .(as.numeric(tas + temp_offset), as.numeric(prec * prec_offset))]
  
  return(as.vector(c(Year_clim = foc_year, 
                     as.matrix(clim_tab_sim[,c("tas","prec","rad", "vpd")]))))
  
}
extract_clim_func <- cmpfun(extract_clim_func)


# wrapper fun
wrapper_fun <- function(x){
  
  x <- as.vector(as.data.frame(x))
  sim_point <- get_point_id(as.data.frame(x[c("Lon", "Lat")]))
  
  simdata_new <- simdata[simulationID == as.numeric(x$ID)
                         & !is.na(Temp) & !is.na(Precip)]
  
  
  # get offset
  temp_offset <- as.numeric(x["temp_offset"])
  prec_offset <- as.numeric(x["prec_offset"])
  
  
  # load the timeseries
  focal_ts <- as.vector(unlist(
    timeseries_focal_df[timeseries_focal_df[, "uniqueID"] == as.numeric(x["uniqueID"]) & 
                          timeseries_focal_df[, "simulationID"] == as.numeric(x["ID"]), ][-c(1,2)]))
  
  simdata_new <- cbind(simdata_new, Year_clim = na.omit(focal_ts))
  
  unique_focal_ts <- unique(na.omit(focal_ts))
  unique_focal_ts_list <- split(unique_focal_ts, seq(length(unique_focal_ts)))  
  
  clim_list <- lapply(unique_focal_ts, FUN = extract_clim_func,
                      clim_tab = clim_tab, sim_point = sim_point,
                      temp_offset = temp_offset, prec_offset = prec_offset)
  
  clims <- t(rbindlist(list(clim_list)))
  
  colnames(clims) <- c("Year_clim", paste0("tas_", 1:365), paste0("prec_", 1:365),
                       paste0("rad_", 1:365), paste0("vpd_", 1:365))
  sim_df <- simdata_new %>% left_join(as.data.frame(clims), by = c("Year_clim"))
  
  return(sim_df)
  
}
wrapper_fun <- cmpfun(wrapper_fun)


### 1. load metadata --------------------------------------------------------------------------
simulation_db <- DBI::dbConnect(
  RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.sqlite"))


tables_con <- dbListTables(simulation_db)
unique(tables_con)
metadata <- dbReadTable(simulation_db, "metadata_climate_scenario")


### 2. order metadata so we can training_db bunch of climate data from spark -------------------
# for instance order with target scenario in order to not switch between the scenario folders 
# and load all grid cells of 1000 simulations at the time
metadata <- metadata %>% arrange(scenario)
metadata <- metadata %>% filter(!is.na(scenario)) %>% data.table()


# load time series for each simulation in bulk, from here we can extract in the loop what we need
timeseries_dfs <- tables_con[grepl("timeseries", tables_con)]

for(i in timeseries_dfs){
  
  tab <- dbReadTable(simulation_db, paste(i)) %>% data.table()
  assign(paste(i), tab)
  
}


### 3. loop over the simulations by groups of scenarios --------------------------------------
scenarios <- metadata %>% dplyr::select(scenario) %>% distinct() %>% unlist() %>% as.vector()
scenarios

for(s in scenarios){
    
    print(s)
    if(is.na(s)){next}
    
    # subset the metadata by scenario
    meta_sub_scen <- metadata %>% fsubset(scenario == s)
    
    # get the coords
    coords <- meta_sub_scen[, c("Lon", "Lat")] %>%
      distinct() %>% data.table()
    
    my_point_id <- apply(coords, 1, get_point_id)
    
    users <- as.vector(unlist(unique(meta_sub_scen[, "uniqueID"])))
    
    ### 4. read spark climate data for each simulation according the scenario and 
    sc <- sparklyr::spark_connect(master = "local", config = config)
    spark_tbl_handle <- spark_read_parquet(sc, name="climate",
                                           path=paste0("/mnt/dss/data/climate_database/spark_db_v3/", s, "/"),
                                           memory = FALSE)
    
    clim_tab <- spark_tbl_handle %>%
      filter(point_id %in% my_point_id) %>% 
      select(., c(point_id, year, month, day, tas, prec, rad, vpd)) %>%
      collect %>%
      setDT()
    
    
    # clear the sparklyr session
    sc %>% spark_session() %>% sparklyr::invoke("catalog") %>% 
      sparklyr::invoke("dropTempView", "spark_tbl_handle")
    sc %>% spark_session() %>% sparklyr::invoke("catalog") %>%
      sparklyr::invoke("clearCache")
    sparklyr::spark_disconnect(sc)
    sparklyr::spark_disconnect_all()
    sc <- NULL
    spark_tbl_handle <- NULL
    rm(sc)
    gc()
    
    
    users_time <- system.time({
      
      for(u in users){
        
        if(!(u %in% c(1031, 1037:1047))){next}
        print(u)
        #if(s == "MPI-M-MPI-ESM-LR_rcp_4_5" & u < 1021){next}
        
        # separate the meta data of the user
        meta_sub_scen_user <- meta_sub_scen %>% fsubset(uniqueID == u)
        
        # load the table with the new timeseries
        timeseries_focal_df <- as.data.frame(get(paste0("user_", u, "_new_timeseries")))
        
        # split to chunks if too large
        chnks <-  split(c(1:nrow(meta_sub_scen_user)),
                        ceiling(seq_along(c(1:nrow(meta_sub_scen_user)))/3000))
        
        # loop over chunks
        for(chnk in 1:length(chnks)){
          
          # print progress  
          print(paste0(chnk," / ", length(chnks)))
          
          # take the metadata of the chunk
          meta_chnk <- meta_sub_scen_user[chnks[[chnk]][1]:chnks[[chnk]][length(chnks[[chnk]])], ]
          
          # convert it to a lst so we can mc lapply over it
          meta_chnk_list <- split(meta_chnk, seq(nrow(meta_chnk)))  
          
          # get all ids
          chnk_ids_sub <- as.vector(unlist(meta_chnk[,"ID"]))
          
          # load the simulation data of all those IDs
          simulation_db <- DBI::dbConnect(
            RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.sqlite"))
          
          
          tables_con <- dbListTables(simulation_db)
          sim_tables <- tables_con[grepl(paste0("simulation"), tables_con)]
          
          simdata <- simulation_db %>%
            tbl(sim_tables[grepl(u, sim_tables)]) %>%
            filter(simulationID %in% chnk_ids_sub) %>%
            collect() %>%
            arrange(simulationID)
          
          simdata <- setDT(simdata)
          
          dbDisconnect(simulation_db)
          
          simdata_new_list <- parallel::mclapply(meta_chnk_list, FUN = wrapper_fun,
                                                 mc.cores = 36, mc.cleanup = TRUE)
          
          gc(verbose = F)
          
          sim_df_user <- rbindlist(simdata_new_list)
          
          # write to DB
          training_db <- DBI::dbConnect(
            RSQLite::SQLite(),
            paste0(path, "/harmonized_db_", s, "_v1.sqlite"))
          
          
          if(chnk == 1){
            dbWriteTable(conn = training_db, name = paste0("harmonized_simulations_", u, "_", s),
                         value = sim_df_user, overwrite = T)
          }else{
            dbWriteTable(conn = training_db, name = paste0("harmonized_simulations_", u, "_", s),
                         value = sim_df_user, append = T)
          } 
          
          
          dbDisconnect(training_db)
          
          sim_df_user <- NULL
          simdata_new_list <- NULL
          simdata <- NULL
          meta_chnk_list <- NULL
          meta_chnk <- NULL
          training_db <- NULL
          simulation_db <- NULL
          rm(simdata, sim_df_user, meta_chnk, meta_chnk_list,
             simdata_new_list, training_db, simulation_db,
             chnk_ids_sub, tables_con, sim_tables)
          
          tmp_dir <- tempdir()
          temp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
          file.remove(temp_files)
          
          for (i in 1:10){gc(reset = T)} 
          
        } # close chnks
        
        
        rm(timeseries_focal_df, meta_sub_scen_user)
        for (i in 1:10){gc(reset = T)} 
        
      } # close user
    })
    
    print(paste0("all user time: ",  users_time[3]))
    
    
    clim_tab <- NULL
    rm(clim_tab, meta_sub_scen)
    for (i in 1:10){gc(reset = T)} 
    
    
}

