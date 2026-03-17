############### 2. RESULTS PT1 #####################

# Creates figures 1, 2, and 3
# Summary results of the extent of seal-fishery interactions
# Kernel utilisation distributions of seal locations and seal-fishery overlap locations
# Gear types involved in seal overlap 

################ CONTENTS ################ 

# A. LIBRARIES, WD, DATASETS
# B. FIGURE 1
# C. OVERLAP SUMMARY
# D. FIGURE 2 
# E. FIGURE 3


################ A. LIBRARIES, WD, DATASETS ################ 

library(tidyverse)
library(lubridate) 
library(hms) 
library(data.table) 
library(sf) 
library(ggplot2) 
library(ggmap) 
library(ggspatial)
library(rnaturalearthdata) 
library(rnaturalearth)
library(adehabitatHR) 
library(patchwork) 
library(move)
library(gfwr)

# Set WD to R Scripts and Data folder
setwd('')

# Country outlines for maps
world <- ne_countries(scale = "large", returnclass = "sf", continent = c("Europe"))

WGS1984 <- CRS("EPSG:4326")
GFW_TOKEN <- ''

# Load datasets
full_interpolated_overlap <- read_csv('full_interpolated_overlap.csv')
missing_time_rows <- grepl("^\\d{4}-\\d{2}-\\d{2}$", full_interpolated_overlap$DateTime)
full_interpolated_overlap$DateTime[missing_time_rows] <- paste0(full_interpolated_overlap$DateTime[missing_time_rows], " 00:00:00")
full_interpolated_overlap <- full_interpolated_overlap %>% 
  mutate(Seal_ID = as.factor(Seal_ID),
         vessel_overlap = as.factor(vessel_overlap),
         gillnet_overlap = as.factor(gillnet_overlap),
         trawl_overlap = as.factor(trawl_overlap),
         foraging = as.factor(foraging),
         DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:00", tz = "UTC"))

overlap_interp_1km <- read_csv('overlap_interp_1km.csv')
missing_time_rows <- grepl("^\\d{4}-\\d{2}-\\d{2}$", overlap_interp_1km$DateTime)
overlap_interp_1km$DateTime[missing_time_rows] <- paste0(overlap_interp_1km$DateTime[missing_time_rows], " 00:00:00")
rm(missing_time_rows)
overlap_interp_1km <- overlap_interp_1km %>% 
  mutate(gear_type = as.factor(gear_type), 
         Seal_ID = as.factor(Seal_ID),
         DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:00", tz = "UTC"))

metadata <- read_csv('metadata_all.csv')

################# B. FIGURE 1 ####################### 

# Map of seal tracks, coloured by individual, with tagging locations added
(track_plot <- ggplot(data = world) + 
  geom_path(data = full_interpolated_overlap, aes(x = Longitude, y = Latitude, colour = Seal_ID), 
            alpha = 0.8, linewidth = 0.4) +
  geom_sf(color = "black", fill = "gray83") +
  annotate("point", y = 50.235262, x = 1.543531, colour = "black", 
           size = 3, alpha = 0.7) + # Baie de Somme
  annotate("point", y = 49.432388, x = 0.137526, colour = "black", fill = 'white', 
           size = 3, alpha = 0.7) + #Seine Estuary
  coord_sf(xlim = c(min(full_interpolated_overlap$Longitude) - 1.8, max(full_interpolated_overlap$Longitude) + 0.5), 
           ylim = c(min(full_interpolated_overlap$Latitude) - 0.5, max(full_interpolated_overlap$Latitude) + 0.5), 
           expand = TRUE) +
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank(), 
        axis.title = element_blank(), axis.text = element_text(size = 11))+
  annotation_scale(data = full_interpolated_overlap, width_hint = 0.4))

ggsave("all_track_plot.png", plot = track_plot)

############## C. OVERLAP SUMMARY ######################

# Number of interpolated locations per individual  
GPS_per_seal <- full_interpolated_overlap %>%
  group_by(Seal_ID) %>%
  summarise(
    n_GPS = n(), hours_tracked = (n() * 14) / 60 # *14 because 14 mins between interpolated locations
  )
mean(GPS_per_seal$hours_tracked) 
sd(GPS_per_seal$hours_tracked) 

# Hours seals were tracked at sea
(nrow(full_interpolated_overlap)*14)/60 

# Hours at sea co-occurring with fishing activity 
nrow(filter(full_interpolated_overlap, vessel_overlap == "Y")) * 14/60 
# Percentage of track co-occurring with fishing activity 
(nrow(filter(full_interpolated_overlap, vessel_overlap == "Y")) / nrow(full_interpolated_overlap))*100
# For each individual, the percentage of track overlapping with fishing gear
seal_ratios <- full_interpolated_overlap %>%
  group_by(Seal_ID) %>%
  summarise(
    overlap_percentage = (sum(vessel_overlap == "Y") / n()) * 100
  )
summary(seal_ratios$overlap_percentage) 
sd(seal_ratios$overlap_percentage) 

# Overlap with just gillnets
nrow(filter(full_interpolated_overlap, gillnet_overlap == "Y"))*14/60 # Number of hours
(nrow(filter(full_interpolated_overlap, gillnet_overlap == "Y"))/nrow(full_interpolated_overlap))*100 # percent of time at sea
# For each individual, the percentage of track overlapping with fishing gear
gillnet_ratios <- full_interpolated_overlap %>%
  group_by(Seal_ID) %>%
  summarise(
    overlap_percentage = (sum(gillnet_overlap == "Y") / n()) * 100
  )
summary(gillnet_ratios$overlap_percentage) 
sd(gillnet_ratios$overlap_percentage) 

# Number of individual vessels overlapping with 
length(unique(overlap_interp_1km$vessel_id)) 
vessel_number_gear <- overlap_interp_1km %>% 
  group_by(gear_type) %>% 
  summarise(vessels = length(unique(vessel_id)))%>%
  mutate(gear_type = recode(gear_type, DREDGE_FISHING = "Dredge", FIXED_GEAR = "Fishing",
                            POTS_AND_TRAPS = "Pots", SET_GILLNETS = "Set gillnets", 
                            TRAWLERS = "Trawl", SEINERS = "Seine", OTHER_SEINES = "Seine",
                            OTHER_PURSE_SEINES = "Seine", FISHING = "Fishing")) %>% 
  aggregate(vessels ~ gear_type, sum)


# Length and number of vessel visits
# Adding foraging to overlap dataframe, then removing duplicate DateTime values of seal at the same vessel 
overlap_filtered <- overlap_interp_1km %>% 
  left_join(full_interpolated_overlap, by = join_by(DateTime, Seal_ID)) %>% 
  dplyr::select(-c(Longitude.y, Latitude.y, vessel_overlap, gillnet_overlap, hourly)) %>% 
  rename(Longitude = Longitude.x, Latitude = Latitude.x) %>% 
  arrange(Seal_ID, DateTime) %>% 
  distinct(Seal_ID, DateTime, vessel_name, .keep_all = TRUE)
# For each vessel interaction, how many consecutive periods of overlap (less than 30 mins between points)
overlap_filtered_2 <- overlap_filtered %>% 
  arrange(Seal_ID, DateTime) %>% 
  group_by(Seal_ID, vessel_name) %>%
  mutate(within_30_mins = ifelse(difftime(DateTime, lag(DateTime), units = "mins") <= 15, "Y", NA)) %>% 
  arrange(.by_group = TRUE, DateTime) %>%
  ungroup() %>% 
  mutate(y_flag = within_30_mins == "Y", 
         run_id = with(rle(y_flag), rep(seq_along(lengths), lengths)), 
         within_30_mins = ifelse(is.na(within_30_mins), "N", within_30_mins), 
         run_id = ifelse(within_30_mins == "N" & lead(within_30_mins) == "Y",
                         lead(run_id), run_id)) 
# Visits per month 
metadata <- metadata %>% 
  mutate(months = Duration/30)
visits_seals <- visits_seals %>% 
  left_join(metadata, by = join_by(Seal_ID)) %>% 
  mutate(visits_month = no_visits / months) %>% 
  dplyr::select(Seal_ID, no_visits, visits_month)
mean(visits_seals$visits_month) 
sd(visits_seals$visits_month) 
# Length of each individual visit
avg_durations <- overlap_filtered_2 %>%
  group_by(Seal_ID, run_id) %>%
  summarise(run_duration = sum(as.numeric(difftime(max(DateTime), min(DateTime), units = "mins")) + 14)) #14 mins added to account for time of final point overlap (missed when taking time difference)
mean(avg_durations$run_duration) 
sd(avg_durations$run_duration)  
# Variation between individuals
avg_y_durations <- overlap_filtered_2 %>%
  group_by(Seal_ID, run_id) %>%
  summarise(run_duration = sum(as.numeric(difftime(max(DateTime), min(DateTime), units = "mins")) + 14)) %>%
  ungroup %>% 
  group_by(Seal_ID) %>%
  summarise(avg_y_duration = mean(run_duration))
mean(avg_y_durations$avg_y_duration) 
sd(avg_y_durations$avg_y_duration) 

# Number of overlap points that are in day/night
nrow(filter(full_interpolated_overlap, sun_elevation > -6, vessel_overlap == "Y"))/
  nrow(filter(full_interpolated_overlap, sun_elevation > -6))*100 # percent during day 
nrow(filter(full_interpolated_overlap, sun_elevation < -6, vessel_overlap == "Y"))/
  nrow(filter(full_interpolated_overlap, sun_elevation < -6))*100 # percent at night 

# Foraging behaviour stats 
# Percentage of track spent exhibiting foraging behaviour 
(nrow(filter(full_interpolated_overlap, foraging == "Y"))/nrow(full_interpolated_overlap))*100 
# Percentage by individual
foraging_percentage <- full_interpolated_overlap %>%
  group_by(Seal_ID) %>%
  summarise(
    foraging_percentage = (sum(foraging == "Y") / n()) * 100
  )
summary(foraging_percentage$foraging_percentage) 
sd(foraging_percentage$foraging_percentage) 

# Foraging during co-occurrence with fishing activity 
# Percentage of whole track
(nrow(filter(full_interpolated_overlap, foraging == "Y", vessel_overlap == "Y")) / 
    nrow(filter(full_interpolated_overlap, foraging == "Y")))*100 #overlapping when foraging
(nrow(filter(full_interpolated_overlap, foraging == "Y", vessel_overlap == "Y")) / 
    nrow(filter(full_interpolated_overlap, vessel_overlap == "Y")))*100 #foraging when overlap
# Variation between individuals
foraging_o_percentage <- full_interpolated_overlap %>% #foraging taking place during overlap
  group_by(Seal_ID) %>%
  summarise(
    foraging_overlap_percent = (sum(foraging == "Y" & vessel_overlap == "Y") / sum(vessel_overlap == "Y") * 100
    ))
summary(foraging_o_percentage$foraging_overlap_percent) 
sd(foraging_o_percentage$foraging_overlap_percent) 

foraging_o_percentage <- full_interpolated_overlap %>% #overlap taking place during foraging
  group_by(Seal_ID) %>%
  summarise(
    foraging_overlap_percent = (sum(foraging == "Y" & vessel_overlap == "Y") / sum(foraging == "Y") * 100
    ))
summary(foraging_o_percentage$foraging_overlap_percent) 
sd(foraging_o_percentage$foraging_overlap_percent) 



################### D. FIGURE 2 ##########################
# A five-panelled figure of (A) map of seal tracks and haulout locations, (B) heat map of fishing activity, 
# (C) spatial distribution of seal at-sea locations, (D) spatial distribution of seal-fishery co-occurrence, 
# (E) distribution of individual seals' proportion of time overlapping with fishing activity

# Map of tracks and haulout locations (panel A) 
haulout_rounded <- haulout_all %>% 
  mutate(across(c('LON', 'LAT'), round, 1)) %>% 
  group_by(LON, LAT) %>% 
  summarise(number_haulouts = n()) %>% 
  na.omit() %>% 
  filter(number_haulouts > 1) %>% 
  mutate(percent = (number_haulouts / nrow(haulout_all))*100)
GPS_haulout <- ggplot(data = world) + 
  geom_path(data = GPS_all, aes(x = Longitude, y = Latitude), 
            alpha = 0.8, linewidth = 0.2, colour = "#98c2fd") +
  geom_sf(color = "black", fill = "gray83") +
  geom_point(data = haulout_rounded, 
             aes(x = LON, y = LAT, 
                 colour = number_haulouts, 
                 size = number_haulouts)) +
  scale_colour_gradient(low = "blue", high = "red", 
                        breaks = c(10, 100, 250, 500, 1000),
                        labels = c(10, 100, 250, 500, 1000),
                        guide = "legend") +
  scale_size(range = c(0.5, 5), 
             breaks = c(10, 100, 250, 500, 1000),
             labels = c(10, 100, 250, 500, 1000),
             guide = "legend") +
  coord_sf(xlim = c(min(GPS_all$Longitude) - 1.8, max(GPS_all$Longitude)) + 0.5, 
           ylim = c(min(GPS_all$Latitude), max(GPS_all$Latitude)) + 0.5, 
           expand = TRUE) +
  theme_bw() +
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(), legend.position = "inside", 
        legend.position.inside = c(0.75,0.84), legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size=10), legend.text = element_text(size=10),
        axis.text = element_text(size = 11)) +
  labs(colour = "Haul-out events", 
       size = "Haul-out events")

ggsave("GPS_haulout.png", GPS_haulout, width = 3.5, height = 5.5)

# Map of fishing activity (panel B) 
gridded_activity <- read.csv("") # downloaded from GFW

gridded_plot <- ggplot()+
  geom_tile(data = gridded_activity, 
            aes(x = Lon, y = Lat, fill = Apparent.Fishing.Hours)) +
  scale_fill_gradient2(low = "white", mid = "red", high = "red", midpoint = 6000,
                       breaks = c(2000, 4000, 6000, 8000),
                       labels = c(2000, 4000, 6000, 8000),
                       guide = "legend")+
  geom_sf(data = world, fill = "grey80") +
  coord_sf(xlim = c(-5.45, 4.28), 
           ylim = c(49.15, 58.15))+
  theme_bw() +
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(), legend.position = "inside", 
        legend.position.inside = c(0.75,0.84), legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size=10), legend.text = element_text(size=10),
        axis.text = element_text(size = 11))+
  labs(fill = "Hours / 120" ~km^2)
ggsave("GFW_effort.png", gridded_plot, width = 3.5, height = 5.5)


# Kernel UD of at-sea locations (panel C) 
gps_utm <- full_interpolated_overlap
coordinates(gps_utm) <- ~Longitude+Latitude 
proj4string(gps_utm) <- CRS("+proj=longlat +datum=WGS84")
gps_utm <- spTransform(gps_utm, CRS("+proj=utm +zone=29 +datum=WGS84"))
kud1 <- kernelUD(gps_utm, h = 15000, grid = 500)
hr95.1 <- getverticeshr(kud1, percent = 95) %>% 
  st_as_sf() %>% 
  st_transform(st_crs(world))
hr50.1 <- getverticeshr(kud1, percent = 50) %>% 
  st_as_sf() %>% 
  st_transform(st_crs(world))
location_kud <- ggplot() +
  geom_sf(data = hr95.1, aes(fill = "95%"), alpha = 0.4, color = "black") +
  geom_sf(data = hr50.1, aes(fill = "50%"), alpha = 0.4, color = "black") +
  geom_sf(data = world, fill = "grey80") +
  coord_sf(xlim = c(min(full_interpolated_overlap$Longitude) - 1.88, max(full_interpolated_overlap$Longitude)) + 0.5, 
           ylim = c(min(full_interpolated_overlap$Latitude), max(full_interpolated_overlap$Latitude)) + 0.5, 
           expand = TRUE) +
  theme(panel.grid = element_blank(), plot.title = element_text(size = 9))+
  scale_fill_manual(name = 'Seal Kernel UD',
                    breaks=c('95%', '50%'),
                    values=c('95%'='lightgreen', '50%'='darkgreen'))+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "inside", 
        legend.position.inside = c(0.8,0.844), legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'), 
        legend.title = element_text(size=10), legend.text = element_text(size=10),
        axis.text = element_text(size = 11))
location_kud

ggsave("location_kud.png", location_kud, width = 3.5, height = 5.5)

# Kernel UD of seal-fishery overlap location (panel D)
overlap_points <- filter(full_interpolated_overlap, vessel_overlap == "Y")
coordinates(overlap_points) <- ~Longitude+Latitude 
proj4string(overlap_points) <- CRS("+proj=longlat +datum=WGS84")
overlap_points_utm <- spTransform(overlap_points, CRS("+proj=utm +zone=29 +datum=WGS84"))
kud2 <- kernelUD(overlap_points_utm, h = 15000, grid = 500) #h is bandwidth
hr95.2 <- getverticeshr(kud3, percent = 95) %>% 
  st_as_sf() %>% 
  st_transform(st_crs(world))
hr50.2 <- getverticeshr(kud3, percent = 50) %>% 
  st_as_sf() %>% 
  st_transform(st_crs(world))
overlap_kud <- ggplot() +
  geom_sf(data = hr95.2, aes(fill = "95%"), alpha = 0.4, color = "black") +
  geom_sf(data = hr50.2, aes(fill = "50%"), alpha = 0.4, color = "black") +
  geom_sf(data = world, fill = "grey80") +
  coord_sf(xlim = c(min(full_interpolated_overlap$Longitude) - 1.88, max(full_interpolated_overlap$Longitude)) + 0.5, 
           ylim = c(min(full_interpolated_overlap$Latitude), max(full_interpolated_overlap$Latitude)) + 0.5, 
           expand = TRUE) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.75,0.844),
        legend.key.size = unit(0.4, 'cm'), legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm'), legend.title = element_text(size=10),
        legend.text = element_text(size=10), axis.text = element_text(size = 11))+
  scale_fill_manual(name = 'Overlap Kernel UD',
                    breaks=c('95%', '50%'),
                    values=c('95%'='lightblue', '50%'='darkblue'))
overlap_kud
ggsave("overlap_kud.png", overlap_kud, width = 3.5, height = 5.5)

# Individual variation in overlap (panel E)
seal_ratios <- full_interpolated_overlap %>%
  group_by(Seal_ID) %>%
  summarise(
    overlap_percentage = (sum(vessel_overlap == "Y") / n()) * 100
  )
overlap_variation_plot <- ggplot(seal_ratios, aes(y = overlap_percentage))+
  geom_boxplot(width = 0.7) +
  xlim(-1,1)+
  labs( y = "%")+
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 12), 
        panel.grid = element_blank(), axis.text.y = element_text(size = 12))
ggsave("overlap_variation_plot.png", overlap_variation_plot, width = 2, height = 5.5)


##################### E. FIGURE 3 ###############################
# A three-panelled figure of barcharts of gear types involved in seal overlap. (A) percentage of overlap time 
# made up by different gear types, (B) Percentage of individual overlapping vessels associated with each gear type,
# (C) number of fishing vessels of each gear type active within the seal tracking region and period

length(unique(overlap_interp_1km$vessel_id)) 
unique_seals <- unique(overlap_interp_1km$Seal_ID)

# Panel A
gear_ratios <- overlap_interp_1km %>%
  group_by(gear_type) %>%
  summarize(overlap_points = n(), .groups = "drop") %>%
  mutate(ratio = (overlap_points / nrow(overlap_interp_1km)*100)) %>%
  dplyr::select(gear_type, ratio) %>% 
  mutate(gear_type = recode(gear_type, DREDGE_FISHING = "Dredge", FIXED_GEAR = "Unknown",
                            POTS_AND_TRAPS = "Pots", SET_GILLNETS = "Set gillnets", 
                            TRAWLERS = "Trawl", SEINERS = "Seine", OTHER_SEINES = "Seine",
                            OTHER_PURSE_SEINES = "Seine", FISHING = "Unknown")) %>% 
  aggregate(ratio ~ gear_type, sum)

# Panel B
vessel_number_gear <- overlap_interp_1km %>% 
  group_by(gear_type) %>% 
  summarise(vessels = length(unique(vessel_id))/158*100)%>%
  mutate(gear_type = recode(gear_type, DREDGE_FISHING = "Dredge", FIXED_GEAR = "Unknown",
                            POTS_AND_TRAPS = "Pots", SET_GILLNETS = "Set gillnets", 
                            TRAWLERS = "Trawl", SEINERS = "Seine", OTHER_SEINES = "Seine",
                            OTHER_PURSE_SEINES = "Seine", FISHING = "Unknown")) %>% 
  aggregate(vessels ~ gear_type, sum)

# Panel C
effort_all <- list()
for(seal in unique_seals){
  #defining aeqd 
  track <- filter(full_interpolated_overlap, Seal_ID == seal)
  Grey_mv <-  move(x = track$Longitude, y = track$Latitude,
                   time = track$DateTime, proj = WGS1984,
                   data = track)
  Grey_mv <- spTransform(Grey_mv, center = TRUE)
  aeqd <- Grey_mv@proj4string
  #track point projection
  track_points <- track %>% 
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>% 
    st_transform(crs = aeqd)
  
  # Create polygon of buffer around track (5km)
  track_buffer <- st_union(st_buffer(track_points, dist = 5000, endCapStyle = "ROUND")) %>% 
    st_transform(crs = 4326) %>%
    st_sf
  
  start_date <- date(min(track$DateTime))
  end_date <- date(max(track$DateTime))
  
  #within a 5km buffer of the tracks, what vessels there on a monthly scale?
  effort <- get_raster(spatial_resolution = 'LOW',
                       temporal_resolution = 'MONTHLY',
                       group_by = 'VESSEL_ID', 
                       start_date = start_date,
                       end_date = end_date,
                       region = track_buffer,
                       region_source = 'USER_SHAPEFILE',
                       key = GFW_TOKEN)
  effort_all[[seal]] <- effort
}

effort_all_df <- bind_rows(effort_all) %>%
  distinct(Lat, Lon, `Time Range`, `Vessel ID`, .keep_all = TRUE) %>% 
  filter(`Gear Type` != 'CARGO', `Gear Type` != 'PASSENGER') %>% 
  mutate(`Gear Type` = recode(`Gear Type`, DREDGE_FISHING = "Dredge", FIXED_GEAR = "Unknown",
                              POTS_AND_TRAPS = "Pots", SET_GILLNETS = "Set gillnets", 
                              TRAWLERS = "Trawl", SEINERS = "Seine", OTHER_SEINES = "Seine",
                              OTHER_PURSE_SEINES = "Seine", FISHING = "Unknown", 
                              TUNA_PURSE_SEINES = "Seine", INCONCLUSIVE = "Unknown",
                              POLE_AND_LINE = "Unknown", SET_LONGLINES = "Longline", 
                              OTHER = "Unknown", PURSE_SEINES = 'Seine')) %>% 
  group_by(`Vessel ID`, `Gear Type`) %>% 
  summarise(fishing_hours = sum(`Apparent Fishing Hours`)) 

effort_all_df <- effort_all_df%>% 
  rename(gear_type = `Gear Type`, vessel_ID = `Vessel ID`) %>% 
  mutate(gear_type = as.factor(gear_type))
rm(effort_all)

length(unique(effort_all_df$vessel_ID)) 
effort_summ <- effort_all_df %>% 
  ungroup() %>% 
  group_by(gear_type) %>% 
  summarise(gear_proportion = n()/nrow(effort_all_df)*100)

palette <- c("#b4ddd4", "#42908c", "#10eddc", "#214d4e", "#66de78", "#3e9539", "#c5df72")
names(palette) <- levels(factor(c(levels(vessel_number_gear$gear_type), levels(gear_ratios$gear_type),
                                  levels(effort_summ$gear_type))))

# Panel A- percentage of overlap time made up by different gear types
vessel_gear_time <- ggplot(gear_ratios, aes(x = reorder(gear_type, -ratio), y = ratio,
                                            fill = gear_type))+
  geom_bar(stat = "identity")+
  theme_bw()+
  labs(x = "Gear Type", y = "%")+
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(name = "gear_type", values = palette)+
  ylim(0, 65)+
  guides(fill = FALSE)

# Panel B- gear type of all individual vessels involved in seal overlap
vessel_gear_IDs <- ggplot(vessel_number_gear, aes(x = reorder(gear_type, -vessels), y = vessels,
                                              fill = gear_type))+
  geom_bar(stat = "identity")+
  theme_bw()+
  labs(x = "Gear Type", y = "%")+
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(name = "gear_type", values = palette)+
  ylim(0, 65)+
  guides(fill = 'none')

# Panel C- all vessels in the study area 
gear_all_plot <- ggplot(effort_summ, aes(x = reorder(gear_type, -gear_proportion), y = gear_proportion,
                                         fill = gear_type))+
  geom_bar(stat = "identity")+
  labs(x = "Gear Type", y = "%")+
  theme_bw()+
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(name = "gear_type", values = palette)+
  ylim(0, 65)+
  guides(fill = FALSE)

# Patchwork the three above together
patch_2 <- (vessel_gear_time / vessel_gear_IDs / gear_all_plot) +
  plot_annotation(tag_levels = 'A')
ggsave(plot = patch_2, filename = "patch_2.png", width = 7.5, height = 10)



