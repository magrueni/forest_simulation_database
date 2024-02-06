### explore dataset ---

library(DBI)
library(tidyverse)


# set puth to dataset
path <- "/your/path/"

# open connection
simulation_db <- DBI::dbConnect(RSQLite::SQLite(),
                                paste0(path, "/forest_simulation_db_v1.sqlite"))

# read tables
tables_con <- dbListTables(simulation_db)
tables_con



### explore metadata ---
metadata <- dbReadTable(simulation_db, "metadata_climate_scenario")

# how many simulations
nrow(metadata)

# number of simulations from X model
models <- unique(metadata$Model)

metadata_model <- metadata %>% 
  filter(Model == models[1])

# explore subset
nrow(metadata_model)





### explore simulation data ---

# loop over all simulation tables and count simulations
# ! careful, this will take a while !
total_sims <- c()

for(i in 1:length(unique(metadata$uniqueID))){
  
  sim_dat <- dbReadTable(simulation_db, paste0("simulationdata_user_", unique(metadata$uniqueID)[i], "_states"))
  total_sims <- c(total_sims, length(unique(sim_dat$simulationID)))
  rm(sim_dat)
  
}

print(sum(total_sims))



# look at a simulation in detail
uniqueID <- 1003 # change uniqueID to look at different simulations
sim_dat <- dbReadTable(simulation_db, paste0("simulationdata_user_1003_states")) 

# select a simulation ID
sim_id <- sample(unique(sim_dat$simulationID), 1) # just to randomly look at one

# select from all simulations of the unique ID
sim_dat %>% 
  filter(simulationID == sim_id) %>% 
  View()



### explore vegetation states ---
states <- dbReadTable(simulation_db, "uniqueStates")

# check species
states <- states %>% 
  rowwise() %>% 
  mutate(species = strsplit(svd_state, "_")[[1]][1],
         LAI = strsplit(svd_state, "_")[[1]][2],
         CanopyHeight = strsplit(svd_state, "_")[[1]][3]) # here the lower bound of the class

states %>% group_by(species) %>% 
  summarize(count = n()) %>% 
  ggplot(., aes(y = count)) +
  geom_bar() +
  ggtitle("Histogram of Sample Data") +
  xlab("Values") +
  ylab("Frequency")



dbDisconnect(simulation_db)





