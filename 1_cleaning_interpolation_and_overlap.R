############# 1. DATA CLEANING, INTERPOLATION AND FISHERY OVERLAP ##############

# Loads seal tracking, dive, and haulout data 
# Interpolates GPS data, removes fixes during haulout behaviour
# Adds columns needed for analysis: bathymetry, day/night, substrate, dive complexity
# Identifies seal locations co-occurring with fishing activity


################ CONTENTS ################ 

# A. LIBRARIES AND WD
# B. LOAD + CLEAN DATASETS
# C. GPS INTERPOLATION
# D. OVERLAP FUNCTION
# E. FORAGING BEHAVIOUR AND OVERLAP
# F. SEDIMENT TYPE
# G. BATHYMETRY
# H. SOLAR POSITION
# I. COMPLEXITY 


################ A. LIBRARIES AND WD ################ 

library(tidyverse)
library(lubridate) 
library(hms) 
library(data.table) 
library(move) 
library(sf) 
library(gfwr) 
#library(pak), pak::pak("EMODnet/emodnet.wfs")
library(emodnet.wfs) 
library(suntools) 
library(marmap)

# Set WD to R Scripts and Data folder  
setwd('')

# Set coordinate reference system, required for 'move' functions
WGS1984 <- CRS("EPSG:4326")

# Save GFW API Token, https://globalfishingwatch.org/our-apis/documentation#quick-start
GFW_TOKEN <- ""


################ B. LOAD + CLEAN DATASETS ################ 

# Read in GPS data, removing fixes with implausible speeds, individuals with tracking for less than a month
GPS_all <- read_csv2('') %>% # Add in GPS CSV file name
  rename(DateTime = D_DATE, Latitude = LAT, Longitude = LON, Seal_ID = Seal) %>% 
  dplyr::select(Seal_ID, DateTime, Latitude, Longitude, V_MASK) %>%
  mutate(
    Seal_ID = as.factor(Seal_ID), V_MASK = as.factor(V_MASK),
    DateTime = as.POSIXct(DateTime, format = "%d/%m/%Y %H:%M", tz = "UTC")
  )%>%
  filter(
    V_MASK == c("0"), # Removing fixes with implausible speeds
    grepl('G', Seal_ID) # Only including grey seals
  ) %>% 
  group_by(Seal_ID) %>% 
  mutate(timediff = difftime(max(DateTime), min(DateTime), units = 'days')) %>% 
  filter(timediff > 30) %>% 
  ungroup() %>% droplevels() %>% 
  dplyr::select(Seal_ID, DateTime, Latitude, Longitude)

# Read in haulout data
haulout_all <- read_csv2('') %>% # Add in haulout CSV file name
  rename(Seal_ID = Seal, Latitude = lat, Longitude = lon) %>%
  mutate(
    Seal_ID = as.factor(Seal_ID),
    S_DATE = as.POSIXct(S_DATE, format = "%d/%m/%Y %H:%M", tz = "UTC"),
    E_DATE = as.POSIXct(E_DATE, format = "%d/%m/%Y %H:%M", tz = "UTC")
  ) %>%
  filter(grepl('G', Seal_ID), Seal_ID %in% GPS_all$Seal_ID) %>% 
  droplevels()


# Read in dive data
dives_all <- read_csv2('') %>% # Add in dive CSV file name
  rename(Seal_ID = Seal, start_lat = lat, start_lon = lon) %>% 
  dplyr::select(-c(REF, PCA_MAX_DESC, PCA_MAX_BTM, PCA_MAX_ASC, PCA_MEAN_ASC, 
                   PCA_MEAN_BTM, PCA_MEAN_DESC, end_lat, end_lon)) %>% 
  drop_na(start_lat, start_lon)%>% 
  mutate(
    start_lat = as.numeric(start_lat),
    start_lon = as.numeric(start_lon),
    Seal_ID = as.factor(Seal_ID),
    DE_DATE = as.POSIXct(DE_DATE, format = "%d/%m/%Y %H:%M", tz = "UTC"),
    MAX_DEP = as.numeric(-MAX_DEP),
    foraging = as.factor(ifelse(PCA_DESC == 1 | PCA_BTM == 1 | PCA_ASC == 1, "Y", "N")),
    PCA_TOTAL = PCA_DESC + PCA_BTM + PCA_ASC, 
    DS_DATE = as.POSIXct(format((DE_DATE - DIVE_DUR), format = '%Y-%m-%d %H:%M:00'), tz = 'UTC'), 
    end_lat = NA, end_lon = NA
  ) %>% 
  filter(grepl('G', Seal_ID), Seal_ID %in% GPS_all$Seal_ID) %>% 
  droplevels()


#trimming GPS and dive data so they match up in tracking period
dive_end_times <- dives_all %>%
  group_by(Seal_ID) %>%
  summarise(max_dive_end = max(DE_DATE, na.rm = TRUE), .groups = "drop")
GPS_all <- GPS_all %>%
  inner_join(dive_end_times, by = "Seal_ID") %>%  # Join to match seals with dive end times
  filter(DateTime <= max_dive_end) %>%            # Keep rows only up to dive end time
  dplyr::select(-max_dive_end)
GPS_end_times <- GPS_all %>%
  group_by(Seal_ID) %>%
  summarise(max_GPS_end = max(DateTime, na.rm = TRUE), .groups = "drop")
dives_all <- dives_all %>%
  inner_join(GPS_end_times, by = "Seal_ID") %>%  # Join to match seals with GPS end times
  filter(DE_DATE <= max_GPS_end) %>%            # Keep rows only up to GPS end time
  dplyr::select(-max_GPS_end)


# Load metadata file
metadata <- read_csv('metadata_all.csv')


############## C. GPS INTERPOLATION #############
# Splitting track when > 90 min gap between fixes to avoid interpolating in long gaps

# Function to split the dataframe based on time difference (time_mins), filtering out single point segments
split_on_time_diff <- function(df, time_mins) {
  #Calculate time difference between points
  df <- df %>%
    mutate(time_diff = c(NA, difftime(DateTime[-1], DateTime[-n()], units = "mins")))
  #Identify indices where time_diff is greater than specified interval
  split_indices <- which(df$time_diff > time_mins)
  #Create split points for slicing the data
  split_points <- c(1, split_indices, nrow(df) + 1)
  
  #Split the dataframe, add segment ID, filter out single-point segments
  split_dfs <- lapply(seq(length(split_points) - 1), function(i) {
    segment_df <- df[split_points[i]:(split_points[i + 1] - 1), ]
    
    # Check if segment has more than one point before adding it to the list
    if (nrow(segment_df) > 1) {
      segment_df <- segment_df %>%
        mutate(segment_id = i)
      return(segment_df)
    } else {
      return(NULL)  
    }
  })
  # Remove <1-point segments
  split_dfs <- split_dfs[!sapply(split_dfs, is.null)]
  return(split_dfs)
}

# Apply the split function to each dataframe in the list, splitting when gap of >90 mins
split_df <- GPS_all %>% 
  split(f = GPS_all$Seal_ID) %>% 
  lapply(split_on_time_diff, time_mins = 90) %>% 
  bind_rows()

# Checking the time difference (TD) between the GPS fixes
add_time_diff <- function(df) {
  mutate(df, time_diff = c(NA, difftime(DateTime[-1], DateTime[-n()], units = "mins"))) 
}
timediff_df <- add_time_diff(GPS_all)
timediff_df$time_diff[timediff_df$Seal_ID != lag(timediff_df$Seal_ID)] <- NA # don't calculate between individuals
summary(timediff_df$time_diff) # median time difference between fixes = 14 mins

# Function to interpolate tracks
interpolate.GPS <- function(df, mins){
  unique_segments <- unique(df$segment_id)
  interpolated_segments <- list() #initialise list
  
  for (segment in unique_segments) {
    # Subset the track for the current segment
    segment_df <- filter(df, segment_id == segment)
    # Convert the track segment to a Move object
    move_object <- move(x = segment_df$Longitude, 
                        y = segment_df$Latitude, 
                        time = segment_df$DateTime, 
                        proj = WGS1984,
                        data = segment_df)
    
    # Apply the interpolation function within the segment
    interpolated_segment <- interpolateTime(move_object, 
                                            time = as.difftime(mins, units = "mins"), 
                                            spaceMethod = "greatcircle")
    
    # Convert the interpolated result to a dataframe and store in list
    interpolated_segments[[segment]] <- as.data.frame(interpolated_segment)
  }
  interpolated_df <- bind_rows(interpolated_segments)
  
  interpolated_df <- interpolated_df %>% 
    dplyr::select(Seal_ID, timestamps, coords.x1, coords.x2, sensor, segment_id) %>% 
    rename(Longitude = coords.x1, Latitude = coords.x2, DateTime = timestamps, Interpolated = sensor) %>%
    mutate(Interpolated = ifelse(Interpolated == "interpolateTime", "Interpolated", "Real"),
           DateTime = as.POSIXct(DateTime, tz = "UTC"))
}

# Applying interpolation at median frequency of time difference between fixes (14 mins here)
full_interpolated_df <- split_df %>% 
  #arrange(Seal_ID, DateTime) %>% 
  split(f = split_df$Seal_ID) %>% 
  lapply(interpolate.GPS, mins = 14) %>% 
  bind_rows()


# Removing any locations during haulout period as only interested in time at sea
full_interpolated_no_haulout <- full_interpolated_df %>%
  mutate(row = row_number()) %>%
  anti_join(haulout_all, by = join_by(Seal_ID, DateTime > S_DATE, DateTime < E_DATE)) %>%
  dplyr::select(-c(Interpolated, segment_id, row))


############### D. OVERLAP FUNCTION ###################

# Function finding overlap for each seal location with active and static fishing activity
# Outputs a dataframe of seal locations that overlap with fishing activity and the nature of that activity (gear type, vessel ID, etc)
# Includes a spatial and temporal buffer to specify when defining seal-fishery overlap
get.point.based.overlap <- function(track, sp_buffer_m, tp_buffer_hr){
  # Create point layer of track, project into Azimuthal Equi-Distance Projection (aeqd) so map units in meters
  Grey_mv <-  move(x = track$Longitude, y = track$Latitude,
                   time = track$DateTime, proj = WGS1984,
                   animal = as.factor(track$Seal_ID), data = track,
                   removeDuplicatedTimestamps = TRUE)
  Grey_mv <- spTransform(Grey_mv, center = TRUE)
  aeqd <- Grey_mv@proj4string
  #track point projection
  track_points <- track %>% 
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>% 
    st_transform(crs = aeqd) %>% 
    mutate(Longitude = track$Longitude, Latitude = track$Latitude,
           hourly = as.POSIXct(format(DateTime, format = '%Y-%m-%d %H:00:00'), tz = "UTC"))
  
  # Create polygon of buffer around track (spatial buffer + 2km)
  track_buffer <- st_union(st_buffer(track_points, dist = sp_buffer_m + 2000, endCapStyle = "ROUND")) %>% 
    st_transform(crs = 4326) %>%
    st_sf
  
  # Create vector of days tracked, find min and max
  start_date <- as.character(date(min(track$DateTime, na.rm = TRUE)))
  end_date <- as.character(date(max(track$DateTime, na.rm = TRUE)))
  
  # Download GFW apparent fishing effort during those days within shapefile
  effort <- get_raster(spatial_resolution = 'HIGH',
                       temporal_resolution = 'HOURLY',
                       group_by = 'VESSEL_ID', 
                       start_date = start_date,
                       end_date = end_date,
                       region = track_buffer,
                       region_source = 'USER_SHAPEFILE',
                       key = GFW_TOKEN)
  # Changing effort to sf object
  effort_points <- effort %>% 
    st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE) %>% 
    st_transform(crs = aeqd) %>% 
    mutate(`Gear Type` = as.factor(`Gear Type`)) %>% 
    rename(gear_type = `Gear Type`, time_vessel = `Time Range`, vessel_id = `Vessel ID`, 
           flag = Flag, vessel_name = `Vessel Name`, gear_type = `Gear Type`, 
           apparent_fishing_hours = `Apparent Fishing Hours`)
  
  #looking for distance overlap on hourly scale 
  all_joined_data <- lapply(unique(track_points$hourly), function(time) {
    #track points of +/- temporal buffer for all gear types, filter later
    track_point_time <- track_points %>% 
      filter(DateTime > (time - tp_buffer_hr*3600) & DateTime < (time + tp_buffer_hr*3600)) %>% 
      dplyr::select(Seal_ID, DateTime, Longitude, Latitude, hourly)
    effort_point_time <- effort_points %>% 
      filter(time_vessel == time) %>% 
      rename(lon_vessel = Lon, lat_vessel = Lat) %>% 
      dplyr::select(time_vessel, lon_vessel, lat_vessel, vessel_id, flag, vessel_name, 
                    gear_type, MMSI, apparent_fishing_hours)
    
    #join track points with effort points if within spatial buffer distance
    joined_data <- st_join(track_point_time, effort_point_time, 
                           join = st_is_within_distance, dist = sp_buffer_m, left = FALSE)
    if (nrow(joined_data) > 0) return(joined_data) else return(NULL)
  })
  # Combine all results into a single data frame, remove duplicate rows
  final_joined_data <- bind_rows(all_joined_data) %>% 
    st_drop_geometry() %>% 
    distinct()
  
  #filter for active gears so only points overlapping on the same hour included
  final_joined_data <- final_joined_data[-which(final_joined_data$gear_type != 'SET_GILLNETS' & 
                                                  final_joined_data$gear_type != 'POTS_AND_TRAPS' & 
                                                  final_joined_data$hourly != final_joined_data$time_vessel),]
  
  return(final_joined_data)
}

# Apply overlap function to interpolated tracks- here with 1km overlap, and a 24 hr buffer for static gear
overlap_interp_1km <- full_interpolated_no_haulout %>% 
  split(f = full_interpolated_no_haulout$Seal_ID) %>% 
  lapply(get.point.based.overlap, sp_buffer_m = 1000, tp_buffer_hr = 24) %>% 
  bind_rows()

# Matching up overlap with track- create column for fishery overlap (Y/N)
unique_seals <- unique(haulout_all$Seal_ID)

full_interpolated_overlap <- full_interpolated_no_haulout %>% 
  mutate(vessel_overlap = NA)

for(seal in unique_seals){ 
  track <- filter(full_interpolated_overlap, Seal_ID == seal)
  overlap <- filter(overlap_interp_1km, Seal_ID == seal)
  full_interpolated_overlap[(full_interpolated_overlap$Seal_ID == seal),] <- track %>%
    rowwise() %>%  
    mutate(
      vessel_overlap = ifelse(any(DateTime == overlap$DateTime),
                              "Y", "N")
    ) %>%
    ungroup()
}

# Matching up overlap with track- create column for gillnet overlap, and for trawl overlap
full_interpolated_overlap <- full_interpolated_no_haulout %>% 
  mutate(gillnet_overlap = NA, trawl_overlap = NA)
for(seal in unique_seals){
  track <- filter(full_interpolated_overlap, Seal_ID == seal)
  overlap <- filter(overlap_interp_1km, Seal_ID == seal, gear_type == "SET_GILLNETS")
  full_interpolated_overlap[(full_interpolated_overlap$Seal_ID == seal),] <- track %>%
    rowwise() %>%  
    mutate(
      gillnet_overlap = ifelse(any(DateTime == overlap$DateTime),
                               "Y", "N")
    ) %>%
    ungroup()
}
for(seal in unique_seals){
  track <- filter(full_interpolated_overlap, Seal_ID == seal)
  overlap <- filter(overlap_interp_1km, Seal_ID == seal, gear_type == "TRAWLERS")
  full_interpolated_overlap[(full_interpolated_overlap$Seal_ID == seal),] <- track %>%
    rowwise() %>%  # Apply operation row by row
    mutate(
      trawl_overlap = ifelse(any(DateTime == overlap$DateTime),
                             "Y", "N")
    ) %>%
    ungroup()
}


############### E. FORAGING BEHAVIOUR AND OVERLAP ####################
# Matching up GPS and dive data with periods of foraging and vessel overlap
# Adding columns to denote this

# Matching up GPS locations with foraging behaviour (i.e. if a prey capture dive occurs in a 14 minute window of that location)
full_interpolated_overlap$foraging <- NA
for(seal in unique_seals){
  seal_track <- filter(full_interpolated_overlap, Seal_ID == seal)
  seal_dives <- filter(dives_all, Seal_ID == seal, foraging == "Y")
  full_interpolated_overlap[(full_interpolated_overlap$Seal_ID == seal),] <- seal_track %>%
    rowwise() %>%  
    mutate(
      foraging = ifelse(nrow(filter(seal_dives, DE_DATE >= (DateTime - 7*60), 
                                    DE_DATE < (DateTime + 7*60))) > 0,
                        "Y", "N")
    )
}


# Adding whether near fishing activity to dives (Y/N) 
# (i.e. if a location overlapping with fishing activity occurs in a 14 minute window of a dive)
dives_all$vessel_overlap <- NA
for(seal in unique_seals){
  dives <- filter(dives_all, Seal_ID == seal)
  overlap <- filter(overlap_interp_1km, Seal_ID == seal)
  dives_all[(dives_all$Seal_ID == seal),] <- dives %>%
    rowwise() %>%  
    mutate(
      vessel_overlap = ifelse(any(abs(difftime(DS_DATE, overlap$DateTime, units = "mins")) <= 7),
                              "Y", "N")
    ) %>%
    ungroup()
}

# Adding whether near gillnet activity to dives (Y/N), and whether near trawl activity
dives_all <- dives_all %>% 
  mutate(gillnet_overlap = NA, trawl_overlap = NA)
for(seal in unique_seals){
  dives <- filter(dives_all, Seal_ID == seal)
  overlap <- filter(overlap_interp_1km, Seal_ID == seal, gear_type == "SET_GILLNETS")
  dives_all[(dives_all$Seal_ID == seal),] <- dives %>%
    rowwise() %>% 
    mutate(
      gillnet_overlap = ifelse(any(abs(difftime(DS_DATE, overlap$DateTime, units = "mins")) <= 7),
                               "Y", "N")
    ) %>%
    ungroup()
}
for(seal in unique_seals){
  dives <- filter(dives_all, Seal_ID == seal)
  overlap <- filter(overlap_interp_1km, Seal_ID == seal, gear_type == "TRAWLERS")
  dives_all[(dives_all$Seal_ID == seal),] <- dives %>%
    rowwise() %>%  
    mutate(
      trawl_overlap = ifelse(any(abs(difftime(DS_DATE, overlap$DateTime, units = "mins")) <= 7),
                             "Y", "N")
    ) %>%
    ungroup()
}


# Removing dives during periods without tracking data to increase certainty of dive location
# All dive points need to have an interpolated location within 14 mins (median difference between )
dives_filtered <- list()
for(seal in unique_seals){
  seal_dives <- filter(dives_all, Seal_ID == seal)
  seal_track <- filter(full_interpolated_overlap, Seal_ID == seal)
  dives_filtered[[seal]] <- seal_dives %>% 
    rowwise() %>% 
    mutate(include = ifelse(any(abs(difftime(DE_DATE, seal_track$DateTime, units = "mins")) <= 14),
                            "Y", "N")) %>% 
    ungroup() %>% 
    filter(include == "Y") %>% 
    dplyr::select(-include)
}


# Remove duplicate locations in proximity to both gillnets and trawls (only 7 dives)
dives_d1 <- dives_all %>% 
  mutate(gillnet_overlap = ifelse(gillnet_overlap == "Y" &
                                  trawl_overlap == "Y", 'N', gillnet_overlap))
dives_d2 <- dives_all %>% 
  filter(gillnet_overlap == 'Y' & trawl_overlap == 'Y') %>% 
  mutate(trawl_overlap = 'N')
dives_all <- rbind(dives_d1, dives_d2)  

# Creating dataset of just foraging dives, add column denoting type of fishery overlap, remove overlap dives near to other gear types
dives_fg_ov <- filter(dives_all, foraging == "Y") %>% 
  mutate(overlap_gear = ifelse(gillnet_overlap == "Y", "Gillnet", "None"),
         overlap_gear = ifelse(trawl_overlap == "Y", "Trawl", overlap_gear),
         overlap_gear = ifelse(vessel_overlap == "Y" & gillnet_overlap == "N" &
                                 trawl_overlap == "N", NA, overlap_gear),
         overlap_gear = as.factor(overlap_gear)) %>% 
  drop_na(overlap_gear) 


############### F. ADDING SEDIMENT TYPE ############### 

# Adding column for sediment type of the location of dives
services <- emodnet_wfs()
seabed_wfs_client <- emodnet_init_wfs_client(
  service = "seabed_habitats_general_datasets_and_products"
)
habitats_directive_layer_names <- "eusm2025_eunis2019_full"

# Creating bounding box of extent of tracks
track_points <- full_interpolated_overlap %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326,
           remove = FALSE)

bbox_vals <- st_bbox(track_points)

bbox_string <- paste(
  bbox_vals["xmin"],
  bbox_vals["ymin"],
  bbox_vals["xmax"],
  bbox_vals["ymax"],
  "EPSG:4326",
  sep = ","
)

# Request seabed habitat map 
habitats_directive_layers <- emodnet_get_layers(
  wfs = seabed_wfs_client,
  layers = habitats_directive_layer_names,
  bbox = bbox_string,
  simplify = TRUE,
  outputFormat = "JSON"
)
habitats_directive_layers_1 <- bind_rows(habitats_directive_layers) %>% 
  dplyr::select(substrate, geometry) %>% 
  st_transform(st_crs(track_points)) # align coordinate systems

# Adding sediment type to dive dataset 
dives_fg_ov <- dives_fg_ov %>% 
  st_as_sf(coords = c("start_lon", "start_lat"), crs = 4326, remove = FALSE)
dives_fg_ov_sub <- st_join(dives_fg_ov, habitats_directive_layers_1) %>% 
  mutate(substrate = as.factor(
    recode(substrate, Sand = "Fine", `Coarse substrate` = "Coarse", `Rock or other hard substrata` = "Rocky", 
           `Muddy sand` = "Fine", `Sandy mud` = "Fine", `Worm reefs` = "Coarse", `Mixed sediment` = "Coarse",
           `[Sabellaria spinulosa] reefs` = "Coarse", Seabed = NA_character_, `[Modiolus modiolus] beds` = "Coarse", 
           `[Mytilus edulis] beds` = "Coarse", Sediment = NA_character_, `Fine mud` = "Fine", 
           `Sandy mud or Muddy sand` = "Fine")
    ))

# Add the substrate type of the closest known location to NA locations
dives_subs <- filter(dives_fg_ov_sub, !is.na(substrate))
nearest_index <- st_nearest_feature(dives_fg_ov_sub, dives_subs)
missing_rows <- which(is.na(dives_fg_ov_sub$substrate))
dives_fg_ov_sub$substrate[missing_rows] <- dives_subs$substrate[nearest_index[missing_rows]]


############## G. BATHYMETRY ################ 

# Add bathymetry to dive data 
bathy <- getNOAA.bathy(lon1 = min(dives_all$start_lon) - 0.1, 
                       lat1 = min(dives_all$start_lat) - 0.1,
                       lon2 = max(dives_all$start_lon) + 0.1, 
                       lat2 = max(dives_all$start_lat) + 0.1,
                       resolution = 1/5)
dives_all <- dives_all %>% 
  mutate(bathymetry = get.depth(bathy, x = dives_all$start_lon, 
                                y = dives_all$start_lat, locator = FALSE)$depth)

#improving matchup using likely depth of sea floor from seal benthic dives
#solving issue with tides impacting depth at a given location and time
dives_all <- dives_all %>% 
  mutate(hourly = as.POSIXct(format(DS_DATE, format = '%Y-%m-%d %H:00:00'), tz = "UTC"))
for(seal in unique_seals){
  seal_dives <- filter(dives_all, Seal_ID == seal)
  
  start <- as.POSIXct(format(min(seal_dives$DS_DATE), format = '%Y-%m-%d %H:00:00'), tz = "UTC")
  end <- as.POSIXct(format(max(seal_dives$DS_DATE), format = '%Y-%m-%d %H:00:00'), tz = "UTC")
  
  time_seq <- data_frame(hourly = seq.POSIXt(start, end, by = "hour")) %>%
    mutate(hourly_seg = row_number())
  seal_dives2 <- left_join(seal_dives, time_seq, by = join_by(hourly))
  
  for(i in unique(seal_dives2$hourly_seg)){
    dives <- filter(seal_dives2, hourly_seg == i, bathymetry < 0)
    
    if(any(dives$bathymetry - dives$MAX_DEP > 0)){
      dives_2 <- filter(dives, bathymetry - MAX_DEP > 0)
      
      mean <- mean(dives_2$bathymetry - dives_2$MAX_DEP)
      
      seal_dives2[seal_dives2$hourly_seg == i,] <- filter(seal_dives2, hourly_seg == i) %>% 
        mutate(bathymetry = bathymetry - mean)
    }
  }
  dives_all[dives_all$Seal_ID == seal,]$bathymetry <- seal_dives2$bathymetry
}

#omit observations with unrealistic bathymetry (more than 5m shallower than max depth, higher than -1.5m)
dives_fg_ov_correct <- dives_fg_ov_sub %>% 
  filter( bathymetry < -1.5, MAX_DEP - bathymetry > -5) %>%  #remove bathymetry values more than 5m shallower than max depth, higher than -1.5m
  st_drop_geometry() #drop geometry column


############### H. SOLAR POSITION ###############
# Add solar position (day/night) to dives

# Coordinates of dives
dive_points <- dives_fg_ov_correct %>% 
  dplyr::select(start_lon, start_lat) %>% 
  as.matrix()
# Find sun elevation at given location and data/time
sun_pos <- solarpos(dive_points, dives_fg_ov_correct$DE_DATE) %>% 
  as.data.frame()
# Add to dataset 
dives_fg_ov_correct <- dives_fg_ov_correct %>% 
  mutate(sun_elevation = ifelse(sun_pos$V2 >= -6, "Day", "Night"))

# Add dive number to avoid issue of repeated times causing issues with temporal correlation structure
dives_fg_ov_correct <- dives_fg_ov_correct %>% 
  mutate(dive_number = row_number())


################# I. COMPLEXITY ################### 
# Calculating vertical dive complexity 

# Change dive data from wide to long format
dives_long <- dives_fg_ov_correct %>% 
  gather(dive_timepoint, depth_m, D1:D12) %>% 
  mutate(dive_timepoint = gsub('D', '', dive_timepoint)) %>% 
  dplyr::select(-c(T1:T12))
dives_long_2 <- dives_fg_ov_correct %>% 
  gather(dive_timepoint, time_percent, T1:T12) %>% 
  mutate(dive_timepoint = gsub('T', '', dive_timepoint))
dives_long$time_percent <- dives_long_2$time_percent

# Adding time of each dive point 
dives_long <- dives_long %>% 
  mutate(time_dive = DS_DATE + DIVE_DUR*(time_percent/100))
# Ordering by dive number
dives_long <- dives_long[order(dives_long$dive_number),]

# Find just the bottom section of the dive depth points, calculate depth difference between each consecutive depth 
dives_btm <- dives_long %>% 
  group_by(dive_number) %>% 
  filter(time_dive > DS_DATE + SECS_DESC, time_dive < (DS_DATE + SECS_DESC + SECS_BTM)) %>% 
  mutate(depth_diff = abs(depth_m - lag(depth_m)))
dives_btm$depth_diff[is.na(dives_btm$depth_diff)] <- 0

btm_sum <- dives_btm %>% 
  group_by(dive_number) %>% 
  summarise(complexity = sum(depth_diff))

dives_fg_ov_correct <- dives_fg_ov_correct %>% 
  left_join(btm_sum, by = join_by(dive_number)) %>% 
  dplyr::select(-c(dive_number)) %>% 
  mutate(complexity = replace_na(complexity, 0))



#Saving datasets
write.csv(GPS_all, 'GPS_all.csv', row.names = FALSE)
write.csv(dives_all, 'dives_all.csv', row.names = FALSE)
write.csv(haulout_all, 'haulout_all.csv', row.names = FALSE)
write.csv(full_interpolated_overlap, 'full_interpolated_overlap.csv', row.names = FALSE)
write.csv(dives_fg_ov_correct, 'dives_fg_ov_correct.csv', row.names = FALSE)




