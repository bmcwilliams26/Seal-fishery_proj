# Seal-fishery_proj
R code from the analysis of Grey Seal tracking data and Global Fishing Watch fishery data. This study was originally carried out with the following scripts in R V4.4.3 (R Core Team, 2025). 

This is the data and code for:


This repository consists of __ scripts used for ________. These scripts include:


The grey seal location data and detailed data on individuals presented in this study are publicly accessible on the SEANOE website (https://doi.org/10.17882/89715, https://doi.org/10.17882/98773). Dive behaviour data is not publicly accessible. Fishery data is available from Global Fishing Watch (https://globalfishingwatch.org/map).



# The variables in these datasets were:
**Seal GPS data:**

Seal_ID: Seal ID based on the species, location and number of previous individuals fitted with tags

D_DATE: Date and time (UTC 0) of the recorded Fastloc GPS location

LAT: Latitude of the Fastloc GPS location (in decimal degrees)

LON: Longitude of the Fastloc GPS location (in decimal degrees)

V_MASK: SMRU filter of the locations. V_Mask=0 is plausible, V_Mask=20 or 50 is rejected due to unrealistic swimming speed (but may still be true when the seal swims in fast tidal currents!), V_Mask=-1 could not be tested due to the lack of previous/next location


