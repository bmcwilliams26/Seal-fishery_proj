################## 4. REPEAT VISITS ####################
# Is the frequency of seal return visits to foraging sites influenced by co-occurrence with fishing activity?

################ CONTENTS ################ 

# A. LIBRARIES, WD, DATASETS
# B. REPEAT VISIT DATASET
# C. GLMMs

################ A. LIBRARIES, WD, DATASETS ################ 

library(tidyverse)
library(lubridate) 
library(hms) 
library(data.table) 
library(sf) 
library(lme4)
library(DHARMa)

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


########### B. REPEAT VISIT DATASET ###############

track_sf <- filter(full_interpolated_overlap, foraging == "Y") %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
foraging_filtered <- list()
for(seal in unique_seals){
  seal_track <- filter(full_interpolated_overlap, Seal_ID == seal)
  seal_track_end <- date(max(seal_track$DateTime)) - 30 #30 days before the end of the tracking data
  
  seal_foraging_filtered <- filter(track_sf, Seal_ID == seal, date(DateTime) < seal_track_end)
  
  foraging_filtered[[seal]] <- seal_foraging_filtered
}
foraging_filtered <- bind_rows(foraging_filtered)%>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

repeat_all <- list()
foraging_filtered$repeat_visit <- NA
for(seal in unique_seals){
  random_seal <- filter(foraging_filtered, Seal_ID == seal)
  track_seal <- filter(track_sf, Seal_ID == seal)
  
  fg_repeat <- st_join(random_seal, track_seal, 
                       join = st_is_within_distance, dist = 1000, left = FALSE) %>% 
    filter(difftime(DateTime.y, DateTime.x, unit = 'hours') > 24, #limit to just return visits at least a day after
           date(DateTime.y) - date(DateTime.x) <= 30)
  random_seal <- random_seal %>% 
    mutate(repeat_visit = ifelse(DateTime %in% fg_repeat$DateTime.x, 'Y', 'N'))
  foraging_filtered[foraging_filtered$Seal_ID == seal,] <- random_seal
  
  repeat_all[[seal]] <- fg_repeat
}
repeat_all <- bind_rows(repeat_all) 

repeat_summ <- repeat_all %>% 
  st_drop_geometry() %>% 
  group_by(Seal_ID.x, DateTime.x) %>% 
  rename(Seal_ID = Seal_ID.x, DateTime = DateTime.x) %>% 
  summarise(n = n_distinct(unique(date(DateTime.y))), .groups = "drop_last")

foraging_filtered2 <- foraging_filtered %>% 
  st_drop_geometry() %>% 
  dplyr::select(Seal_ID, DateTime, Longitude, Latitude, vessel_overlap, repeat_visit) %>% 
  mutate(vessel_overlap = factor(vessel_overlap),
         vessel_overlap = relevel(vessel_overlap, ref = "N"),
         repeat_visit = ifelse(repeat_visit == 'Y', 1, 0))
foraging_filtered2 <- left_join(foraging_filtered2, repeat_summ,
                                by = join_by(Seal_ID, DateTime))


############### C. GLMMs ################

# repeat probability
repeatprob_glmm <- glmer(repeat_visit ~ vessel_overlap + (1 | Seal_ID),
                data = foraging_filtered2,
                family = binomial())
summary(repeatprob_glmm)
(1 - exp(-0.07335)) *100 #7% less for gear-associated
# effect size
exp(fixef(repeatprob_glmm))
exp(confint(repeatprob_glmm, parm = "beta_", method = "Wald"))

nrow(filter(foraging_filtered2, repeat_visit == 1, vessel_overlap == 'Y'))/
  nrow(filter(foraging_filtered2, vessel_overlap == 'Y')) * 100 #% for gear-associated
nrow(filter(foraging_filtered2, repeat_visit == 1, vessel_overlap == 'N'))/
  nrow(filter(foraging_filtered2, vessel_overlap == 'N')) * 100 #% for not-gear-associated

#LRT
model0 <- glmer(repeat_visit ~ (1 | Seal_ID),
                data = foraging_filtered2,
                family = binomial())

anova(model0, repeatprob_glmm, test = "Chisq")


# Repeat frequency 
foraging_filtered3 <- filter(foraging_filtered2, repeat_visit == 1)

model3 <- glmer(n ~ vessel_overlap + (1|Seal_ID),
                data = foraging_filtered3,
                family = poisson())

repeatfreq_tnb <- glmmTMB(n ~ vessel_overlap + (1 | Seal_ID),
                     data = foraging_filtered2,
                     family = truncated_nbinom2())

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  chisq <- sum(rp^2)
  ratio <- chisq / rdf
  p <- pchisq(chisq, df = rdf, lower.tail = FALSE)
  c(chisq = chisq, ratio = ratio, rdf = rdf, p = p)
}
overdisp_fun(model3) # Testing overdispersion for Poisson model
AIC(model3, model_nb, repeatfreq_tnb)  

simulationOutput<-simulateResiduals(fittedModel = repeatfreq_tnb, plot=T)
testDispersion(repeatfreq_tnb) # Testing overdispersion for tnb model
testOutliers(simulationOutput, type=c("bootstrap")) 
preds<-ggpredict(repeatfreq_tnb)
plot(preds)

#LRT
model_tnb0 <- glmmTMB(n ~ (1 | Seal_ID),
                      data = foraging_filtered2,
                      family = truncated_nbinom2())
anova(model_tnb0, repeatfreq_tnb)

summary(repeatfreq_tnb)
1 - exp(-0.10046)