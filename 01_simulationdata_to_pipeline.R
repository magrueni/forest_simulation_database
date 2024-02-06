##############################################################################
### 01 --- Pipeline for simulation data preparation
##############################################################################
### In this script we bring together all simulation data to a SQLite database
### At the same time we obtain vegetation states from the simulation data
### M. Grünig, 06.02.2024
##############################################################################


library(dplyr)
library(data.table)
library(terra)
library(sf)
library(ggplot2)
library(DBI)
library(stringr)
library(compiler)
library(parallel)


# set memory limit
memory.limit()
memory.size(max=20000)


# path to folder
path <- "/your/path/"
path_responses <- paste0(path, "/your/path/responses")


# list all users
users <- list.files(paste0(path_responses))
users <- users[!(users %in% empty_users)]


### functions ---------------------------------------------------------------------------------------------

source("functions/check_height_cols_function.R")
source("functions/check_attribute_cols_function.R")
source("functions/check_species_cols_function.R")
source("functions/reorder_strings.R")
source("functions/kahn_factors_function_v2.R")
source("functions/height_max_mean_function_v3.R")
source("functions/check_species_condition.R")

### -------------------------------------------------------------------------------------------------------


# mean max models for canopy height harmonization
all_coefs_domheightmodel <- read.csv(paste0(path, "/height_model_coefs.csv"))
all_coefs_maxheightmodel <- read.csv(paste0(path, "/maxheight_model_coefs.csv"))


# create SQLite DB
simulation_db <- DBI::dbConnect(RSQLite::SQLite(),
                                paste0(path, "/forest_simulation_db_v1.sqlite"))


### first look at the metadata and see where the simulations are from -----------------------------------------

# empty objects to store results
metadata_all_users <- NULL
uniqueStates <- NULL


# loop over users to collect simulation data
for(u in users){
    
    print(u)
    userID <- gsub("user_", "", u)
    
    # get the files that belong to the user
    files <- list.files(paste0(path_responses, u))
    metafls <- files[which(grepl("meta", files))]
    
    # create empty lists to store results
    simulation_ids <- c()
    meta_all <- NULL
    
    # if no metafile is available, skip. Else loop over metafiles (there can be more than 1)
    if(length(metafls) == 0){next}else{
      
      # loop over metafiles
      for(m in 1:length(metafls)){    
        print(paste0(m, "/", length(metafls))) 
        
        # read metafile in
        metadata <- read.csv(paste0(path_responses, u, "/", metafls[m]), sep = ",", header = T)
        if(nrow(metadata) != nrow(metadata %>% distinct())){cat("ERROR. Duplicated metadata removed!")}
        metadata <- metadata %>% distinct()
        
        # do some tests to make sure the data has the right format. Remove first column as some data has rownumbers as first col
        if(colnames(metadata)[1] == "V1"){metadata <- metadata[,-c(1)]}
        if(!("uniqueID" %in% colnames(metadata))){metadata <- cbind(uniqueID = rep(as.numeric(userID), nrow(metadata)), metadata)}
        if(!colnames(metadata)[1] == "uniqueID"){metadata <- metadata[,-c(1)]}
        
        # bind together the metadata in case there are several files
        meta_all <- rbind(meta_all, metadata)
        
        # collect the simulation IDs in the metafile to check later if the simulation data is represented in the metafiles
        IDs <- metadata$ID
        if(length(IDs) != length(unique(IDs))){break}
        simulation_ids <- c(simulation_ids, IDs)
        
      }
      
      # check whether there are any duplicates in the metafiles
      if(length(simulation_ids) != length(unique(simulation_ids))){print("Duplicated metafile entry!")}
      
      
      # make sure all columns we need are there. If not and the colum is not required to have values, we add it as NAs
      colnms <- c("uniqueID", "ID", "Model", "ModelDOI", "Lon", "Lat", "Country", "WHC", "TextureSand",
                  "TextureSilt", "TextureClay", "SoilDepth", "SoilWaterRating",
                  "AvailableNitrogen",  "FertilityRating", "Climate", "GCM", "RCM", "Management")
      
      
      meta_new <- NULL
      for(c in 1:length(colnms)){
        
        if(colnms[c] %in% colnames(meta_all)){meta_new <- cbind(meta_new, meta_all[, colnms[c]])}else{
          meta_new <- cbind(meta_new, rep(NA, nrow(meta_all)))
        }
      }
      
      # format to numeric
      meta_new <- as.data.frame(meta_new)
      colnames(meta_new) <- colnms
      meta_new <- meta_new %>% mutate(uniqueID = as.numeric(uniqueID),
                                      ID = as.numeric(ID),
                                      Lon = as.numeric(Lon),
                                      Lat = as.numeric(Lat),
                                      WHC = as.numeric(WHC),
                                      TextureSand = as.numeric(TextureSand),
                                      TextureSilt = as.numeric(TextureSilt),
                                      TextureClay = as.numeric(TextureClay),
                                      SoilDepth = as.numeric(SoilDepth),
                                      SoilWaterRating = as.numeric(SoilWaterRating),
                                      AvailableNitrogen = as.numeric(AvailableNitrogen),
                                      FertilityRating = as.numeric(FertilityRating)
                                      
      )
      
      metadata_all_users <- rbind(metadata_all_users, meta_new)
      
    }
    
    ### simulation data -------------------------------------------------------------
    simfls <- files[which(grepl("simulation", files))]
    
    sim_ids <- c()
    all_sims <- NULL
    
    for(m in 1:length(simfls)){
      print(paste0(m, "/", length(simfls)))  
      simdata <- fread(paste0(path_responses, u, "/", simfls[m]), header = T)
      simdata <- check_species_cols(simdata)
      simdata <- check_height_cols(simdata)
      simdata <- check_attribute_cols(simdata)
      if("V1" %in% colnames(simdata)){simdata <- simdata[,-c("V1")]}
      if(!("uniqueID" %in% colnames(simdata))){simdata <- cbind(uniqueID = rep(as.numeric(userID), nrow(simdata)), simdata)}
      if(!("simulationID" %in% colnames(simdata)) & "ID" %in% colnames(simdata)){
        colnames(simdata)[colnames(simdata) == "ID"] <- "simulationID" }
      if(!("simulationID" %in% colnames(simdata))){simdata %>% mutate(simulationID = ID)}
      
      
      if(as.numeric(userID) %in% c(1012, 1031:1047)){
        
        # check if wierd simulations are in there
        df_list <- split(simdata, simdata$simulationID)
        result_list <- unlist(as.vector(mclapply(df_list, check_species_condition, mc.cores = 24)))
        sims_to_remove <- which(result_list)
        if(length(sims_to_remove) > 0){
          #print(paste0(length(sims_to_remove), "/", length(df_list), " simulations removed"))
          final_df_list <- df_list[-sims_to_remove]
          simdata_filtered <- do.call(rbind, final_df_list)
        }else{
          simdata_filtered <- simdata
        }
        
        all_sims <- rbind(all_sims, simdata_filtered)
        
        rm(df_list, result_list, sims_to_remove, final_df_list)
        
      }else{
        
        all_sims <- rbind(all_sims, simdata)
        
      }
      
    }
    
    rm(simdata)
    
    
    # check that all entries are unique
    if(nrow(all_sims) != nrow(all_sims %>% distinct())){cat("ERROR. Duplicated simulation data!")}
    all_sims_distinct <- all_sims %>% distinct()
    
    if(nrow(all_sims) != nrow(all_sims_distinct)){print("Duplicated simulations!")}
    
    # check if all sim IDs are in the metafile
    if(length(which(!unique(all_sims$ID) %in% unique(meta_all$ID))) != 0){cat("simulations do not match metadata")}
    
    # get sim ID
    simulation_identification <- unique(all_sims_distinct$simulationID)
    rm(all_sims)
    
 
    ### convert data to vegetation states --------------------------------------------
    
    # sort out potential NAs in LAI
    all_sims_distinct <- all_sims_distinct[!is.na(all_sims_distinct$LAI),]
    
    # check for mispelling of pinus sylvestris...
    all_sims_distinct <- all_sims_distinct %>% 
      mutate(Species1 = gsub("pisi", "pisy", Species1)) %>% 
      mutate(Species2 = gsub("pisi", "pisy", Species2)) %>% 
      mutate(Species3 = gsub("pisi", "pisy", Species3)) %>% 
      mutate(Species4 = gsub("pisi", "pisy", Species4)) %>% 
      mutate(Species5 = gsub("pisi", "pisy", Species5))
    
    
    
    all_sims_distinct_states <- all_sims_distinct %>% mutate(Proportion1 = ifelse(is.na(Proportion1), 0, as.numeric(Proportion1)),
                                                             Proportion2 = ifelse(is.na(Proportion2), 0, as.numeric(Proportion2)),
                                                             Proportion3 = ifelse(is.na(Proportion3), 0, as.numeric(Proportion3)),
                                                             Proportion4 = ifelse(is.na(Proportion4), 0, as.numeric(Proportion4)),
                                                             Proportion5 = ifelse(is.na(Proportion5), 0, as.numeric(Proportion5))) %>% 
      mutate(sp_state = paste0("")) %>% 
      mutate(sp_state = ifelse(Proportion1 >= 66, paste0(toupper(Species1)), 
                               ifelse(Proportion2 >= 66, paste0(toupper(Species2)), 
                                      ifelse(Proportion3 >= 66, paste0(toupper(Species3)), 
                                             ifelse(Proportion4 >= 66, paste0(toupper(Species4)),
                                                    ifelse(Proportion5 >= 66, paste0(toupper(Species5)), "")))))) %>%
      
      mutate(sp_state2 = ifelse(Proportion1 < 66 & Proportion1 >= 20, paste0(as.character(sp_state), as.character(Species1)), paste0(sp_state))) %>% 
      mutate(sp_state2 = ifelse(Proportion2 < 66 & Proportion2 >= 20, paste0(as.character(sp_state2), as.character(Species2)), paste0(sp_state2))) %>% 
      mutate(sp_state2 = ifelse(Proportion3 < 66 & Proportion3 >= 20, paste0(as.character(sp_state2), as.character(Species3)), paste0(sp_state2))) %>% 
      mutate(sp_state2 = ifelse(Proportion4 < 66 & Proportion4 >= 20, paste0(as.character(sp_state2), as.character(Species4)), paste0(sp_state2))) %>% 
      mutate(sp_state2 = ifelse(Proportion5 < 66 & Proportion5 >= 20, paste0(as.character(sp_state2), as.character(Species5)), paste0(sp_state2)))
    
    #rm(all_sims_distinct)
    gc()
    
    # reorder species combinations alphabetically
    all_sims_distinct_states <- reorder_strings(all_sims_distinct_states)
    all_sims_distinct_states <- all_sims_distinct_states %>% filter(sp_state2 != "")
    
    ### add LAI categories -------------------------------------------------------------   
    all_sims_distinct_states <- all_sims_distinct_states %>% 
      mutate(lai_state = case_when(LAI < 2 ~ "1",
                                   LAI <= 4 ~ "2",
                                   LAI > 4 ~ "3")) 
    
    ### add the Height -----------------------------------------------------------------
    
    # Separate to the three cases
    dominants <- all_sims_distinct_states[!is.na(all_sims_distinct_states$MaxHeight) & str_detect(all_sims_distinct_states$sp_state2,"[[:upper:]]"), ]
    if(nrow(dominants) > 0){dominants <- kahn_factors_function_cmp(dominants)}else{dominants <- NULL}
    
    
    maxima <- data.table(all_sims_distinct_states[!is.na(all_sims_distinct_states$MaxHeight) & !is.na(all_sims_distinct_states$MinHeight) & !str_detect(all_sims_distinct_states$sp_state2,"[[:upper:]]"), ])
    # if(nrow(maxima) > 0){maxima <- maxima %>% mutate(Hdom = apply(., 1, height_max_mean_function))}else{maxima <- NULL}
    if (nrow(maxima) > 0) {
      results <- mclapply(seq_len(nrow(maxima)), function(i) {
        height_max_mean_function(as.data.frame(maxima[i, ]))
      }, mc.cores = 24)
      
      maxima$Hdom <- unlist(as.vector(results))
    } else {
      maxima <- NULL
    }
    
    
    averages <- data.table(all_sims_distinct_states[is.na(all_sims_distinct_states$MaxHeight) & !is.na(all_sims_distinct_states$MeanHeight), ])
    if (nrow(averages) > 0) {
      results <- mclapply(seq_len(nrow(averages)), function(i) {
        height_mean_only_function_cmp(as.data.frame(averages[i, ]))
      }, mc.cores = 24)
      
      averages$Hdom <- unlist(as.vector(results))
    } else {
      averages <- NULL
    }
    
    # combine again, here also rows with no heights are dropped
    all_sims_distinct_states <- bind_rows(dominants, maxima, averages)
    
    # put into states
    all_sims_distinct_states <- all_sims_distinct_states %>% 
      mutate(height_state = case_when(Hdom < 2 ~ "0_2",
                                      Hdom < 4 ~ "2_4",
                                      Hdom < 6 ~ "4_6",
                                      Hdom < 8 ~ "6_8",
                                      Hdom < 10 ~ "8_10",
                                      Hdom < 12 ~ "10_12",
                                      Hdom < 14 ~ "12_14",
                                      Hdom < 16 ~ "14_16",
                                      Hdom < 18 ~ "16_18",
                                      Hdom < 20 ~ "18_20",
                                      Hdom < 22 ~ "20_22",
                                      Hdom < 24 ~ "22_24",
                                      Hdom < 26 ~ "24_26",
                                      Hdom < 28 ~ "26_28",
                                      Hdom < 30 ~ "28_30",
                                      Hdom < 32 ~ "30_32",
                                      Hdom < 34 ~ "32_34",
                                      Hdom < 36 ~ "34_36",
                                      Hdom < 38 ~ "36_38",
                                      Hdom < 40 ~ "38_40",
                                      Hdom < 42 ~ "40_42",
                                      Hdom < 44 ~ "42_44",
                                      Hdom < 46 ~ "44_46",
                                      Hdom < 48 ~ "46_48",
                                      Hdom < 50 ~ "48_50",
                                      Hdom >= 50 ~ "50_")) 
    
    
    ### put together the states --------------------------------------------------------
    all_sims_distinct_states$svd_state <- paste0(all_sims_distinct_states$sp_state2, "_", all_sims_distinct_states$lai_state, "_", all_sims_distinct_states$height_state)

    # sort out the LAI NAs
    all_sims_distinct_states <- all_sims_distinct_states[!is.na(all_sims_distinct_states$LAI)]
    
    # remove columns that are not needed anymore
    all_sims_distinct_states <- all_sims_distinct_states[, -c("sp_state")]
    
    if(length(grep("NA", all_sims_distinct_states$svd_state)) > 0){
      
      cat("NAs detected!")
      cat(all_sims_distinct_states[grep("NA", all_sims_distinct_states$svd_state), "svd_state"], "n\ ")
      
    }
    
    all_sims_distinct_states <- all_sims_distinct_states %>% arrange(simulationID, Year)
  
    uniqueStates_sp <- all_sims_distinct_states %>% dplyr::select(svd_state) %>% dplyr::group_by(svd_state) %>% summarize(count = n())
    uniqueStates <- rbind(uniqueStates, uniqueStates_sp)
    
    dbWriteTable(conn = simulation_db, name = paste0("simulationdata_", u, "_states"), value = all_sims_distinct_states, overwrite = T)
    
  }


# write out   
uniqueStates <- uniqueStates %>% group_by(svd_state) %>% summarize(count = sum(count))
dbWriteTable(conn = simulation_db, name = "uniqueStates", value = uniqueStates, overwrite = T)

dbWriteTable(conn = simulation_db, name = "metadata_all", value = metadata_all_users, overwrite = T)

# disconnect
dbDisconnect(simulation_db)

