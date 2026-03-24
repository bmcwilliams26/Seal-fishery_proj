################## 4. REPEAT VISITS ####################
# Is the frequency of seal return visits to foraging sites influenced by co-occurrence with fishing activity?

################ CONTENTS ################ 

# A. LIBRARIES, WD, DATASETS
# B. FILTERING DATASETS
# C. GEAR-ASSOCIATED RETURN
# D. NON-GEAR-ASSOCIATED RETURN


################ A. LIBRARIES, WD, DATASETS ################ 

library(tidyverse)
library(lubridate) 
library(hms) 
library(data.table) 
library(sf) 

#set working directory 
setwd("")

#metadata
metadata <- read_csv('metadata_all.csv') %>% 
  mutate(Start_tracking = as.POSIXct(Start_tracking, format = '%d/%m/%Y %H:%M', tz = 'UTC'),
         End_tracking = as.POSIXct(End_tracking, format = '%d/%m/%Y %H:%M', tz = 'UTC'))

# Interpolated track
full_interpolated_overlap <- read_csv('full_interpolated_overlap.csv')
#dealing with midnight issue
missing_time_rows <- grepl("^\\d{4}-\\d{2}-\\d{2}$", full_interpolated_overlap$DateTime)
full_interpolated_overlap$DateTime[missing_time_rows] <- paste0(full_interpolated_overlap$DateTime[missing_time_rows], " 00:00:00")
full_interpolated_overlap <- full_interpolated_overlap %>% mutate(Seal_ID = as.factor(Seal_ID),
                                                                  DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:00", tz = "UTC"))
unique_seals = unique(full_interpolated_overlap$Seal_ID)


#################### B. FILTERING DATASETS ################ 
# Creating the dataset for repeat visit analysis following seal-fishery co-occurrence

foraging_overlap <- filter(full_interpolated_overlap, vessel_overlap == "Y", 
                           foraging == "Y") #only foraging locations with overlap
foraging_no_overlap <- filter(full_interpolated_overlap, foraging == "Y", 
                       vessel_overlap == "N") #only foraging locations without overlap

overlap_filtered <- list()
random_filtered <- list()
for(seal in unique_seals){
  seal_track <- filter(full_interpolated_overlap, Seal_ID == seal)
  seal_track_end <- date(max(seal_track$DateTime)) - 30 #30 days before the end of the tracking data
  
  seal_overlap_filtered <- filter(foraging_overlap, Seal_ID == seal, date(DateTime) < seal_track_end)
  seal_random_filtered <- filter(foraging_no_overlap, Seal_ID == seal, date(DateTime) < seal_track_end)
  
  overlap_filtered[[seal]] <- seal_overlap_filtered
  random_filtered[[seal]] <- seal_random_filtered
}

#making tracks into sf objects
overlap_filtered <- bind_rows(overlap_filtered)%>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
random_filtered <- bind_rows(random_filtered)%>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

#sf object of foraging locations
track_sf <- filter(full_interpolated_overlap, foraging == "Y") %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)


#################### C. GEAR-ASSOCIATED RETURN ################ 

repeat_overlap_2 <- list()
for(seal in unique_seals){
  overlap_seal <- filter(overlap_filtered, Seal_ID == seal)
  track_seal <- filter(track_sf, Seal_ID == seal)
  
  fg_repeat <- st_join(overlap_seal, track_seal, #join all points that are within 1km
                       join = st_is_within_distance, dist = 1000, left = FALSE) %>% 
    filter(difftime(DateTime.y, DateTime.x, unit = 'hours') > 24, #limit to just return visits at least a day after 
           date(DateTime.y) - date(DateTime.x) <= 30) #limit to just return within a month
  
  repeat_overlap_2[[seal]] <- fg_repeat
}
repeat_overlap_2 <- bind_rows(repeat_overlap_2) %>% 
  st_drop_geometry()

repeat_overlap_summ <- repeat_overlap_2 %>% 
  group_by(Seal_ID.x, DateTime.x, gillnet_overlap.x, trawl_overlap.x) %>% 
  summarise(n = n_distinct(unique(date(DateTime.y)))) #summarising, number of unique days that return visit

# Chance of a repeat visit
nrow(repeat_overlap_summ) / nrow(overlap_filtered) * 100 
mean(repeat_overlap_summ$n) #average number of repeat visits within a month 
sd(repeat_overlap_summ$n)  


################# D. NON-GEAR-ASSOCIATED RETURN ################# 
# Random subsampling 1000 times of foraging locations without overlap (n = 363)

set.seed(123)

null_percentetc <- data.frame()
for(x in seq(1:1000)){
  locations_fg_sample <- random_filtered %>% 
    sample_n(361) 
  
  repeat_null <- list()
  for(seal in unique_seals){
    random_seal <- filter(locations_fg_sample, Seal_ID == seal)
    track_seal <- filter(track_sf, Seal_ID == seal)
    
    fg_repeat <- st_join(random_seal, track_seal, 
                         join = st_is_within_distance, dist = 1000, left = FALSE) %>% 
      filter(difftime(DateTime.y, DateTime.x, unit = 'hours') > 24, #limit to just return visits at least a day after
             date(DateTime.y) - date(DateTime.x) <= 30)
    
    repeat_null[[seal]] <- fg_repeat
  }
  repeat_null <- bind_rows(repeat_null) 
  
  repeat_null_summ <- repeat_null %>% 
    group_by(Seal_ID.x, DateTime.x) %>% 
    summarise(n = n_distinct(unique(date(DateTime.y))), .groups = "drop_last")
  
  output <- c(nrow(repeat_null_summ)/nrow(locations_fg_sample)*100, mean(repeat_null_summ$n), sd(repeat_null_summ$n))
  null_percentetc <- rbind(null_percentetc, output)
}
colnames(null_percentetc) <- c('percent', 'mean', 'sd')

# Plotting distributions
ggplot(data = null_percentetc, aes(x = percent))+
  geom_histogram(binwidth = 1)
ggplot(data = null_percentetc, aes(x = mean))+
  geom_histogram(binwidth = 1)
ggplot(data = null_percentetc, aes(x = sd))+
  geom_histogram(binwidth = 1)



# Calculating confidence intervals
results <- t.test(null_percentetc$percent)
results$conf.int
results <- t.test(null_percentetc$mean)
results$conf.int

# Empirical two-sided p values 
null_mean <- mean(null_percentetc$percent)
gear_percent <- nrow(repeat_overlap_summ) / nrow(overlap_filtered) * 100 
p_two_sided <- (1 + sum(abs(null_percentetc$percent - null_mean) >= abs(gear_percent - null_mean))) / (nrow(null_percentetc) + 1)
p_two_sided

null_mean_revisit <- mean(null_percentetc$mean)
gear_mean <- mean(repeat_overlap_summ$n)
p_two_sided_mean <- (1 + sum(abs(null_percentetc$mean - null_mean_revisit) >= abs(gear_mean - null_mean_revisit))) / (nrow(null_percentetc) + 1)
p_two_sided_mean

