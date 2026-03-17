# Seal-fishery_proj
R code from the analysis of Grey Seal tracking data and Global Fishing Watch fishery data. This study was originally carried out with the following scripts in R V4.4.3 (R Core Team, 2025). 

This is the data and code for:


This repository consists of __ scripts used for ________. These scripts include:


The grey seal location data and detailed data on individuals presented in this study are publicly accessible on the SEANOE website (https://doi.org/10.17882/89715, https://doi.org/10.17882/98773). Dive behaviour data is not publicly accessible. Fishery data is available from Global Fishing Watch (https://globalfishingwatch.org/map).



# The variables in the seal datasets were:
**Seal GPS data:**

Seal_ID: Seal ID based on the species, location and number of previous individuals fitted with tags

D_DATE: Date and time (UTC 0) of the recorded Fastloc GPS location

LAT: Latitude of the Fastloc GPS location (in decimal degrees)

LON: Longitude of the Fastloc GPS location (in decimal degrees)

V_MASK: SMRU filter of the locations. V_Mask=0 is plausible, V_Mask=20 or 50 is rejected due to unrealistic swimming speed (but may still be true when the seal swims in fast tidal currents!), V_Mask=-1 could not be tested due to the lack of previous/next location


**Dive data:**

Seal_ID: Individual number of the seal fitted with a tag. 

DE_DATE: Date and time of the Dive End

SURF_DUR: Post-dive surface duration (seconds)

DIVE_DUR: Dive duration (seconds)

MAX_DEP: Maximum depth (meters)

D1:Intermediate depth point 1 (metres) in the dive, at time = 10% of the total dive duration
D2:Intermediate depth point 2 (metres) in the dive, at time = 20% of the total dive duration
D3:Intermediate depth point 3 (metres) in the dive, at time = 30% of the total dive duration
D4:Intermediate depth point 4 (metres) in the dive, at time = 40% of the total dive duration
D5:Intermediate depth point 5 (metres) in the dive, at time = 50% of the total dive duration
D6:Intermediate depth point 6 (metres) in the dive, at time = 60% of the total dive duration
D7:Intermediate depth point 7 (metres) in the dive, at time = 70% of the total dive duration
D8:Intermediate depth point 8 (metres) in the dive, at time = 80% of the total dive duration
D9:Intermediate depth point 9 (metres) in the dive, at time = 90% of the total dive duration

LAT:Estimated latitude of the end of the dive (interpolated between previous and next "real" Fastloc GPS location recorded by the tag)

LON:Estimated longitude of the end of the dive (interpolated between previous and next "real" Fastloc GPS location recorded by the tag)


**Haulout data:**

Seal_ID:Telemetry tag number

D_DATE:Date and time (UTC 0) of the recorded Fastloc GPS location

LAT:Latitude of the Fastloc GPS location (in decimal degrees)

LON:Longitude of the Fastloc GPS location (in decimal degrees)

V_MASK:SMRU filter of the locations. V_Mask=0 is plausible, V_Mask=20 is rejected due to unrealistic swimming speed, V_Mask=-1 could not be tested due to the lack of previous/next location
