##############################################################################
### 07 ---  Get the climate database slopes
##############################################################################
### In this script we extract data from the climate database and calculate the
### slope and the averages for all points that we have in the simulation DB
### M. Grünig, 06.02.2024
##############################################################################

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



### spark config ---------------------------------------------------------------

# Set memory allocation for whole local Spark instance
Sys.setenv("SPARK_MEM" = "13g")

# Set driver and executor memory allocations
config <- spark_config()
config$`spark.driver.memory` <- "12G"
config$`sparklyr.shell.driver-memory` <- "12G"
config$`spark.executor.memory` <- "12G"
config$`spark.driver.maxResultSize` <- "12G"
config$`spark.yarn.executor.memoryOverhead` <- "12G"
config$sparklyr.gateway.port = 8892
config$sparklyr.gateway.start.timeout = 180
config$sparklyr.gateway.connect.timeout = 1
config$`sparklyr.cores.local` <- 4
config$`spark.executor.cores` <- 4

options(java.parameters = "-Xmx8048m")



### load reference grid ---
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



### look at the simulations -----------------------------------
simulation_db <- DBI::dbConnect(RSQLite::SQLite(), "/forest_simulation_db_v1.sqlite")
metadata <- dbReadTable(simulation_db, "metadata_all")



coordinates <- metadata[, c("Lon", "Lat")] %>% distinct()

# check which rcp we have to look at
rcps <- c("historical", "rcp_2_6", "rcp_4_5", "rcp_8_5")

# check if the GCM is in there
gcms <- c("ICHEC-EC-EARTH", "MPI-M-MPI-ESM-LR", "NCC-NorESM1-M")


# get the corresponding grid ID for each simulation coordinate
my_point_id <- apply(coordinates, 1, get_point_id)
points_coords_key <- cbind(coordinates, my_point_id)


run_time <- system.time({
  
  # add a counter to know if we are in the first GCM RCP combo or not
  counter <- 0
  
  # loop over gcms and rcps
  for(g in gcms){
    
    for(c in rcps){
      
      counter <- counter + 1  
      
      # connect to spark and open the database of the GCM and RCP combo
      sc <- sparklyr::spark_connect(master = "local", config = config)
      spark_tbl_handle <- spark_read_parquet(sc, name="climate", path=paste0("/mnt/dss/data/climate_database/spark_db_v3/", g, "_", c, "/"), memory = FALSE)
      
      # collect data from spark
      clim_tab <- spark_tbl_handle %>% filter(point_id %in% my_point_id) %>% 
        select(., c(point_id, year, month, day, max_temp, min_temp, prec)) %>% collect() %>% data.table()
      
      # create empty list
      list_slopes_all_points <- list()
      mean_annual_temp <- list()
      annual_prec <- list()
      
      # loop over all points to get temperature and precip data and calculate the slope and the annual average/sum
      r <- 0
      for(p in unique(clim_tab$point_id)){
        r <- r + 1
        
        # filter by point and summarize for years
        clim_tab_p <- clim_tab %>% filter(point_id == p) %>% 
          arrange(year, month, day) %>% 
          group_by(year) %>% summarize(mean_max_temp = mean(max_temp),
                                       mean_min_temp = mean(min_temp),
                                       annual_prec = sum(prec)) %>% 
          arrange(year) %>% 
          mutate(mean_temp = (mean_max_temp + mean_min_temp)/2) %>% 
          mutate(year = as.numeric(year))
        
        # linear models to get slope
        clim_temp_mod <- flm(mean_temp ~ year, data = clim_tab_p)
        clim_prec_mod <- flm(annual_prec ~ year, data = clim_tab_p)
        
        # annual mean temp and annual precip for offset
        clim_temp_avg <- mean(clim_tab_p$mean_temp) 
        clim_prec_avg <- mean(clim_tab_p$annual_prec)
        
        # give everything to an ouput dataframe
        output_df <- data.table(point_id = p, temp_slope = as.numeric(clim_temp_mod[2,1]), prec_slope = as.numeric(clim_prec_mod[2,1]),
                                temp_avg = clim_temp_avg, prec_avg = clim_prec_avg)
        
        
        mean_annual_temp[[r]] <- c(p, t(clim_tab_p[ "mean_temp"]))
        annual_prec[[r]] <- c(p, t(clim_tab_p[ "annual_prec"]))
        
        # add to list
        rownames(output_df) <- r
        list_slopes_all_points[[r]] <- output_df
      }
      
      # bind the list to a dataframe and give names the columns
      all_slopes <- data.table(do.call(rbind, list_slopes_all_points))
      colnames(all_slopes) <- c("point_id", paste0("temp_slope_", g, "_", c), paste0("prec_slope_", g, "_", c), paste0("temp_avg_", g, "_", c), paste0("prec_avg_", g, "_", c))
      
      # add everything to the overall dataframe
      if(counter == 1){df_all <- all_slopes}else{df_all <- cbind(df_all, all_slopes[, 2:5])}
      
      # write the mean annual temp and annual precips
      mean_annual_temp_df <- data.table(do.call(rbind, mean_annual_temp))
      annual_prec_df <- data.table(do.call(rbind, annual_prec))
      
      colnames(mean_annual_temp_df) <- c("my_point_id", as.character(unlist(clim_tab_p[,"year"])))
      mean_annual_temp_df <- left_join(points_coords_key, mean_annual_temp_df, by = c("my_point_id")) %>% rename(., "point_id" = "my_point_id")
      
      colnames(annual_prec_df) <- c("my_point_id", as.character(unlist(clim_tab_p[,"year"])))
      annual_prec_df <- left_join(points_coords_key, annual_prec_df, by = c("my_point_id")) %>% rename(., "point_id" = "my_point_id")
      
      dbWriteTable(conn = simulation_db, name = paste0("annual_mean_temp_", g, "_", c), value = mean_annual_temp_df, overwrite = T)
      dbWriteTable(conn = simulation_db, name = paste0("annual_precip_", g, "_", c), value = annual_prec_df, overwrite = T)
      
      
      
      # clear the sparklyr session
      sc %>% spark_session() %>% sparklyr::invoke("catalog") %>% sparklyr::invoke("dropTempView", "spark_tbl_handle")
      sc %>% spark_session() %>% sparklyr::invoke("catalog") %>% sparklyr::invoke("clearCache")
      sparklyr::spark_disconnect(sc)
      sparklyr::spark_disconnect_all()
      gc()
      
    }
  }
  
})

dbDisconnect(simulation_db)

# join back to metadata in order to make comparable with the simulation data
colnames(points_coords_key) <- c("Lon", "Lat", "point_id") 
head(points_coords_key)

df_all_new <- points_coords_key %>% left_join(df_all, by = c("point_id"))
dbWriteTable(conn = simulation_db, name = paste0("climate_slopes_coordinates"), value = df_all_new, overwrite = T)



