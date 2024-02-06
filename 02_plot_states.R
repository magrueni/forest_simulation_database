##############################################################################
### 02 --- Plotting created vegetation states
##############################################################################
### This script helps to get a graphical overview of the created svd states
### M. Grünig, 06.02.2024
##############################################################################



library(DBI)
library(tidyverse)
library(ggplot2)
library(parallel)
library(doParallel)
library(compiler)
library(data.table)
library(MetBrewer)


path <- "your/path/"


# open database
simulation_db <- DBI::dbConnect(RSQLite::SQLite(),
                                paste0(path, "/forest_simulation_db_v1.sqlite"))

# load states
tables_con <- dbListTables(simulation_db)
unique(tables_con)
states <- dbReadTable(simulation_db, "uniqueStates")
head(states)
dim(states)


# order states by dominant species
states_sp <- states %>% drop_na() %>% 
  mutate(dom_sp = substr(svd_state, 0, 4)) %>% 
  filter(str_detect(dom_sp,"[[:upper:]]")) %>% 
  group_by(dom_sp) %>% 
  summarize(count_states = n(),
            datapoints = sum(count))

ggplot(states_sp) + 
  geom_bar(aes(y = count_states, x = dom_sp, fill = dom_sp), stat="identity") + 
  theme_bw() +
  theme(legend.position = c(x = 20, y = 500))



states_new <- states_sp %>% mutate(dom_sp = as.factor(dom_sp))

ggplot(states_new) + 
  geom_bar(aes(y = datapoints, x = dom_sp),
           stat="identity") +
  scale_y_log10() + 
  theme_bw() +
  theme(legend.position = c(x = 20, y = 500)) +
  scale_fill_discrete() +
  theme(axis.text.y=element_text(size=14),
        axis.title=element_text(size=18),
        axis.text.x=element_text(size=10, angle = -45)) +
  labs(x = "species composition state", y = "number of datapoints")





# order states by LAI
states_lai <- states %>% drop_na() %>% 
  mutate(lai = sapply(str_split(svd_state, "_"), function(x) x[2])) %>% 
  group_by(lai) %>% 
  summarize(count_states = n(),
            datapoints = sum(count))

ggplot(states_lai) + 
  geom_bar(aes(y = count_states, x = lai, fill = lai), stat="identity")

ggplot(states_lai) + 
  geom_bar(aes(y = datapoints, x = lai), fill = rev(met.brewer("Archambault", n = 3)), stat="identity") + 
  theme_bw() +
  theme(legend.position = c(x = 20, y = 500)) +
  scale_fill_discrete() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18)) +
  labs(x = "LAI state", y = "number of datapoints")


# order states by Height
states_height <- states %>% drop_na() %>% 
  mutate(height = sapply(str_split(svd_state, "_"), function(x) x[3])) %>% 
  mutate(height = as.numeric(height)) %>% 
  group_by(height) %>% 
  summarize(count_states = n(),
            datapoints = sum(count))

# plot number of discrete state per height class
ggplot(states_height) + 
  geom_bar(aes(y = count_states, x = height, fill = height), stat="identity")

# plot number of simulation years per height class
ggplot(states_height) + 
  geom_bar(aes(y = datapoints, x = height, fill = height), stat="identity")+ 
  theme_bw() +
  theme(legend.position = c(x = 20, y = 500)) +
  scale_fill_continuous(trans = 'reverse') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18)) +
  labs(x = "dominant height state", y = "number of datapoints")



# all states
all_states <- states %>%
  # filter(count > 1000000) %>% 
  drop_na()

ggplot(all_states) + 
  geom_bar(aes(y = count, x = svd_state), fill = "lightgrey", stat="identity")


