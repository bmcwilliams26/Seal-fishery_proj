############# 3. Linear mixed-effects modelling ###############
# Modelling dive metrics as a response to environmental variables and co-occurrence with fisheries 

################ CONTENTS ################ 

# A. LIBRARIES, WD, DATASETS
# B. DIVE DURATION
# C. TIME AT BOTTOM
# D. SWIM EFFORT
# E. DIVE COMPLEXITY
# F. PLOTS


################ A. LIBRARIES, WD, DATASETS ################ 

library(tidyverse)
library(data.table)
library(lubridate)
library(hms)
library(nlme)
library(sjPlot)

#working directory
setwd("") #insert working directory path 

dives_fg_ov <- read_csv("dives_fg_ov_correct.csv")
# Dealing with issue of midnight values being lost from time columns
missing_time_rows <- grepl("^\\d{4}-\\d{2}-\\d{2}$", dives_fg_ov$DE_DATE)
dives_fg_ov$DE_DATE[missing_time_rows] <- paste0(dives_fg_ov$DE_DATE[missing_time_rows], " 00:00:00")
missing_time_rows <- grepl("^\\d{4}-\\d{2}-\\d{2}$", dives_fg_ov$DS_DATE)
dives_fg_ov$DS_DATE[missing_time_rows] <- paste0(dives_fg_ov$DS_DATE[missing_time_rows], " 00:00:00")
rm(missing_time_rows)

# Changing some from characters to factors
sealDives <- dives_fg_ov %>% 
  mutate(overlap_gear = factor(overlap_gear, levels = c("None", "Gillnet", "Trawl")),
         substrate = factor(substrate, levels = c("Fine", "Coarse", "Rocky")),
         sun_elevation = factor(sun_elevation, levels = c("Day", "Night")),
         DE_DATE = as.POSIXct(DE_DATE, format = "%Y-%m-%d %H:%M:00", tz = "UTC"),
         DS_DATE = as.POSIXct(DS_DATE, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
         Seal_ID = as.factor(Seal_ID),
         gillnet_overlap = as.factor(gillnet_overlap),
         trawl_overlap = as.factor(trawl_overlap)) %>% 
  dplyr::select(Seal_ID, DE_DATE, DIVE_DUR, MAX_DEP, SWIM_EFF_WHOLE, SECS_BTM, 
                PCA_TOTAL, DS_DATE, gillnet_overlap, bathymetry, 
                percent_dep, trawl_overlap, overlap_gear, sun_elevation,
                substrate, wiggliness) %>% 
  rename(sealID = Seal_ID,
         dateTime_start = DS_DATE,
         dateTime_end = DE_DATE,
         maxDepth = MAX_DEP,
         
         diveDur = DIVE_DUR,
         secsBtm = SECS_BTM,
         pca_total = PCA_TOTAL,
         swimEff_total = SWIM_EFF_WHOLE,
         complexity = wiggliness,
         
         overlapGear = overlap_gear,
         percentDepth = percent_dep,
         sunElevation = sun_elevation,
         gillnetOverlap = gillnet_overlap,
         trawlOverlap = trawl_overlap)
#str(sealDives)

#making dive number specific to each individual (needed for autocorrelation function)
sealDives$diveNum <- NA
for(x in unique(sealDives$sealID)){
  seal_data <- filter(sealDives, sealID == x)
  seal_data <- seal_data %>% 
    mutate(diveNum = row_number())
  sealDives[sealDives$sealID == x,] <- seal_data
}
rm(x)
rm(seal_data)

# Add time since last dive, needed for variogram
sealID <- unique(sealDives$sealID)
sealDives$timeSinceStart <- rep(NA, nrow(sealDives))

for(II in 1:length(sealID)){
  indRef <- which(sealDives$sealID == sealID[II])
  sealDives$timeSinceStart[indRef] <- difftime(sealDives$dateTime_start[indRef],
                                               min(sealDives$dateTime_start[indRef]), unit = 'hours')
  rm(indRef)
}
rm(II, sealID)


################ B. DIVE DURATION ####################

#Look at response variable 
plot(sealDives$diveDur)
hist(sealDives$diveDur, breaks = 20)

#Select correlation structure
#base model REML
lme_dur_base <- lme(diveDur ~ overlapGear + substrate + sunElevation + bathymetry,
                    random = ~1|sealID, data = sealDives, method = "REML")
# with temporal correlation structure AR1
lme_dur_ar1 <- lme(diveDur ~ overlapGear + substrate + sunElevation + bathymetry,
                   random = ~1|sealID, data = sealDives, method = "REML",
                   correlation = corAR1(form = ~diveNum|sealID))  
#evaluate models 
AIC(lme_dur_base, lme_dur_ar1)

## lme_dur_base
resid_Model <- resid(lme_dur_base, type='normalized')
plot(resid_Model)

plot(lme_dur_base)
qqnorm(resid_Model)
qqline(resid_Model)

acf(resid_Model) 
plot(Variogram(lme_dur_base, form= ~timeSinceStart|sealID, data = sealDives, resType = 'normalized', robust=TRUE), ylim=c(0,1))

hist(resid_Model, breaks = 40)
boxplot(resid_Model ~ sealDives$sealID, main = "resid ~ sealID",ylim=c(-5,5))

## lme_diveDur_corAR1
resid_Model <- resid(lme_dur_ar1, type='normalized')
plot(resid_Model)

plot(lme_dur_ar1)
qqnorm(resid_Model)
qqline(resid_Model)

#temporal autocorrelation plots
acf(resid_Model) 
plot(Variogram(lme_dur_ar1, form= ~timeSinceStart|sealID, data = sealDives, resType = 'normalized', robust=TRUE), ylim=c(0,1))

hist(resid_Model, breaks = 40)
boxplot(resid_Model ~ sealDives$sealID, main = "resid ~ sealID",ylim=c(-5,5))

# model selection
#testing for the significance of adding each fixed effect in to the model, fitting as maximum likelihood to allow comparison
lme_base <- lme(diveDur ~ overlapGear + substrate + sunElevation + bathymetry,
                random = ~1|sealID, data = sealDives, method = "ML",
                correlation = corAR1(form = ~diveNum|sealID))
lme_nsub <- lme(diveDur ~ overlapGear + sunElevation + bathymetry,
                random = ~1|sealID, data = sealDives, method = "ML",
                correlation = corAR1(form = ~diveNum|sealID))
lme_nov <- lme(diveDur ~ substrate + sunElevation + bathymetry,
               random = ~1|sealID, data = sealDives, method = "ML",
               correlation = corAR1(form = ~diveNum|sealID))
lme_nday <- lme(diveDur ~ overlapGear + substrate + bathymetry,
                random = ~1|sealID, data = sealDives, method = "ML",
                correlation = corAR1(form = ~diveNum|sealID))
lme_nbath <- lme(diveDur ~ overlapGear + substrate + sunElevation ,
                 random = ~1|sealID, data = sealDives, method = "ML",
                 correlation = corAR1(form = ~diveNum|sealID))
tests <- list(
  int   = anova(lme_nov, lme_base),
  add   = anova(lme_nsub, lme_base),
  add2  = anova(lme_nday, lme_base),
  add3 = anova(lme_nbath, lme_base)
)
#table with outputs of likelihood ratio tests
tests2 <- tibble(
  comparison = c("no overlap vs additive",
                 "no substrate vs additive",
                 "no time of day vs additive",
                 "no bathymetry vs additive"),
  L.Ratio = sapply(tests, function(x) x[2, "L.Ratio"]),
  p    = sapply(tests, function(x) x[2, "p-value"])
) # significant addition of all variables
tests2

#final model 
tab1 <- tab_model(lme_dur_ar1, title = 'Dive duration',
                  file = 'divedur.HTML') 



############### C. TIME AT BOTTOM ###############
#correlation structure
#base model 
lme_sbtm_base <- lme(secsBtm ~ overlapGear + substrate + sunElevation + bathymetry,
                     random = ~1|sealID, data = sealDives, method = "REML")
#with sqrt and temporal correlation
lme_sbtm_ar1 <- lme(secsBtm ~ overlapGear + substrate + sunElevation + bathymetry,
                    random = ~1|sealID, data = sealDives, method = "REML",
                    correlation = corAR1(form = ~diveNum | sealID)) 

AIC(lme_sbtm_base, lme_sbtm_ar1)

## lme_sbtm_base ##
resid_Model <- resid(lme_sbtm_base, type='normalized')
plot(resid_Model)

plot(lme_sbtm_base)
qqnorm(resid_Model)
qqline(resid_Model)

acf(resid_Model)
plot(Variogram(lme_sbtm_base, form= ~timeSinceStart|sealID, data = sealDives, resType = 'normalized', robust=TRUE), ylim=c(0,1))

hist(resid_Model, breaks = 40) 
boxplot(resid_Model ~ sealDives$sealID, main = "resid ~ sealID",ylim=c(-5,5))

## lme_sbtm_ar1 ##
resid_Model <- resid(lme_sbtm_ar1, type='normalized')
plot(resid_Model)

plot(lme_sbtm_base)
qqnorm(resid_Model)
qqline(resid_Model)

#temporal autocorrelation plots
acf(resid_Model) 
plot(Variogram(lme_sbtm_ar1, form= ~timeSinceStart|sealID, data = sealDives, resType = 'normalized', robust=TRUE), ylim=c(0,1))

hist(resid_Model, breaks = 40)
boxplot(resid_Model ~ sealDives$sealID, main = "resid ~ sealID",ylim=c(-5,5))

# Model selection
lme_sbtm <- lme(secsBtm ~ overlapGear + substrate + sunElevation + bathymetry,
                random = ~1|sealID, data = sealDives, method = "ML",
                correlation = corAR1(form = ~diveNum|sealID))
lme_sbtm_nsub <- lme(secsBtm ~ overlapGear + sunElevation + bathymetry,
                     random = ~1|sealID, data = sealDives, method = "ML",
                     correlation = corAR1(form = ~diveNum|sealID))
lme_sbtm_nov <- lme(secsBtm ~ substrate + sunElevation + bathymetry,
                    random = ~1|sealID, data = sealDives, method = "ML",
                    correlation = corAR1(form = ~diveNum|sealID))
lme_sbtm_nday <- lme(secsBtm ~ overlapGear + substrate + bathymetry,
                     random = ~1|sealID, data = sealDives, method = "ML",
                     correlation = corAR1(form = ~diveNum|sealID))
lme_sbtm_nbathy <- lme(secsBtm ~ overlapGear + substrate + sunElevation,
                       random = ~1|sealID, data = sealDives, method = "ML",
                       correlation = corAR1(form = ~diveNum|sealID))
tests <- list(
  int   = anova(lme_sbtm_nov,   lme_sbtm),
  add   = anova(lme_sbtm_nsub,   lme_sbtm),
  add2  = anova(lme_sbtm_nday,   lme_sbtm),
  add3 = anova(lme_sbtm_nbathy, lme_sbtm)
)
#table with useful outputs
tests <- tibble(
  comparison = c("no overlap vs additive",
                 "no substrate vs additive",
                 "no time of day vs additive",
                 "no bathymetry vs additive"),
  L.Ratio = sapply(tests, function(x) x[2, "L.Ratio"]),
  p    = sapply(tests, function(x) x[2, "p-value"])
) # all variables are significant additions 
tests

tab2 <- tab_model(lme_sbtm_ar1) 


################ D. SWIM EFFORT #####################

# Removing individual N07 with the dodgy accelerometer 
sealDives_2 <- filter(sealDives, sealID != "N07")

#correlation structure
#base model 
lme_swimeff_base <- lme(swimEff_total ~ overlapGear + substrate + sunElevation + bathymetry,
                        random = ~1|sealID, data = sealDives_2, method = "REML")
#with sqrt and temporal correlation
lme_swimeff_ar1 <- lme(swimEff_total ~ overlapGear + substrate + sunElevation + bathymetry,
                       random = ~1|sealID, data = sealDives_2, method = "REML",
                       correlation = corAR1(form = ~diveNum | sealID)) 

AIC(lme_swimeff_base, lme_swimeff_ar1)

## lme_swimeff_base ##
plot(lme_swimeff_base)
qqnorm(residuals(lme_swimeff_base))
qqline(residuals(lme_swimeff_base))

resid_Model <- resid(lme_swimeff_base, type='normalized')
plot(resid_Model)

acf(resid_Model)  
plot(Variogram(lme_swimeff_base, form= ~timeSinceStart|sealID, data = sealDives_2, resType = 'normalized', robust=TRUE), ylim=c(0,1))

hist(resid_Model, breaks = 40) 
boxplot(resid_Model ~ sealDives_2$sealID, main = "resid ~ sealID",ylim=c(-5,5))

## lme_swimeff_ar1 ##
resid_Model <- resid(lme_swimeff_ar1, type='normalized')
plot(resid_Model)

plot(lme_swimeff_ar1)
qqnorm(resid_Model)
qqline(resid_Model)

#temporal autocorrelation plots
acf(resid_Model) 
plot(Variogram(lme_swimeff_ar1, form= ~timeSinceStart|sealID, data = sealDives_2, resType = 'normalized', robust=TRUE), ylim=c(0,1))

hist(resid_Model, breaks = 40)
boxplot(resid_Model ~ sealDives_2$sealID, main = "resid ~ sealID",ylim=c(-5,5))


# Model selection
lme_swimeff <- lme(swimEff_total ~ overlapGear + substrate + sunElevation + bathymetry,
                   random = ~1|sealID, data = sealDives_2, method = "ML",
                   correlation = corAR1(form = ~diveNum|sealID))
lme_swimeff_nsub <- lme(swimEff_total ~ overlapGear + sunElevation + bathymetry,
                        random = ~1|sealID, data = sealDives_2, method = "ML",
                        correlation = corAR1(form = ~diveNum|sealID))
lme_swimeff_nov <- lme(swimEff_total ~ substrate + sunElevation + bathymetry,
                       random = ~1|sealID, data = sealDives_2, method = "ML",
                       correlation = corAR1(form = ~diveNum|sealID))
lme_swimeff_nday <- lme(swimEff_total ~ overlapGear + substrate + bathymetry,
                        random = ~1|sealID, data = sealDives_2, method = "ML",
                        correlation = corAR1(form = ~diveNum|sealID))
lme_swimeff_nbathy <- lme(swimEff_total ~ overlapGear + substrate + sunElevation,
                          random = ~1|sealID, data = sealDives_2, method = "ML",
                          correlation = corAR1(form = ~diveNum|sealID))
tests <- list(
  int   = anova(lme_swimeff_nov,   lme_swimeff),
  add   = anova(lme_swimeff_nsub,   lme_swimeff),
  add2  = anova(lme_swimeff_nday,   lme_swimeff),
  add3 = anova(lme_swimeff_nbathy, lme_swimeff)
)
#table with useful outputs
tests <- tibble(
  comparison = c("no overlap vs additive",
                 "no substrate vs additive",
                 "no time of day vs additive",
                 "no bathymetry vs additive"),
  L.Ratio = sapply(tests, function(x) x[2, "L.Ratio"]),
  p    = sapply(tests, function(x) x[2, "p-value"])
) # all variables are significant additions 
tests

#final model
tab3 <- tab_model(lme_swimeff_ar1, file = 'swimeff.HTML')


############## E. DIVE COMPLEXITY #############
#correlation structure
#base model 
lme_complex_base <- lme(complexity ~ overlapGear + substrate + sunElevation + bathymetry,
                       random = ~1|sealID, data = sealDives, method = "REML")
#with temporal correlation
lme_complex_ar1 <- lme(complexity ~ overlapGear + substrate + sunElevation + bathymetry,
                      random = ~1|sealID, data = sealDives, method = "REML",
                      correlation = corAR1(form = ~diveNum | sealID)) 

AIC(lme_complex_base, lme_complex_ar1)

## lme_complex_base ##
plot(lme_complex_base)
qqnorm(residuals(lme_complex_base))
qqline(residuals(lme_complex_base))

resid_Model <- resid(lme_complex_base, type='normalized')
plot(resid_Model)

acf(resid_Model) 
plot(Variogram(lme_complex_base, form= ~timeSinceStart|sealID, data = sealDives, resType = 'normalized', robust=TRUE), ylim=c(0,1))

hist(resid_Model, breaks = 40) 
boxplot(resid_Model ~ sealDives_2$sealID, main = "resid ~ sealID",ylim=c(-5,5))

## lme_complex_ar1 ##
resid_Model <- resid(lme_complex_ar1, type='normalized')
plot(resid_Model)

plot(lme_complex_ar1)
qqnorm(resid_Model)
qqline(resid_Model)

#temporal autocorrelation plots
acf(resid_Model) 
plot(Variogram(lme_complex_ar1, form= ~timeSinceStart|sealID, data = sealDives, resType = 'normalized', robust=TRUE), ylim=c(0,1))

hist(resid_Model, breaks = 40)
boxplot(resid_Model ~ sealDives_2$sealID, main = "resid ~ sealID",ylim=c(-5,5))


# Model selection
lme_sbtm <- lme(complexity ~ overlapGear + substrate + sunElevation + bathymetry,
                random = ~1|sealID, data = sealDives, method = "ML",
                correlation = corAR1(form = ~diveNum|sealID))
lme_sbtm_nsub <- lme(complexity ~ overlapGear + sunElevation + bathymetry,
                     random = ~1|sealID, data = sealDives, method = "ML",
                     correlation = corAR1(form = ~diveNum|sealID))
lme_sbtm_nov <- lme(complexity ~ substrate + sunElevation + bathymetry,
                    random = ~1|sealID, data = sealDives, method = "ML",
                    correlation = corAR1(form = ~diveNum|sealID))
lme_sbtm_nday <- lme(complexity ~ overlapGear + substrate + bathymetry,
                     random = ~1|sealID, data = sealDives, method = "ML",
                     correlation = corAR1(form = ~diveNum|sealID))
lme_sbtm_nbathy <- lme(complexity ~ overlapGear + substrate + sunElevation,
                       random = ~1|sealID, data = sealDives, method = "ML",
                       correlation = corAR1(form = ~diveNum|sealID))
tests <- list(
  int   = anova(lme_sbtm_nov,   lme_sbtm),
  add   = anova(lme_sbtm_nsub,   lme_sbtm),
  add2  = anova(lme_sbtm_nday,   lme_sbtm),
  add3 = anova(lme_sbtm_nbathy, lme_sbtm)
)
#table with useful outputs
tests <- tibble(
  comparison = c("no overlap vs additive",
                 "no substrate vs additive",
                 "no time of day vs additive",
                 "no bathymetry vs additive"),
  L.Ratio = sapply(tests, function(x) x[2, "L.Ratio"]),
  p    = sapply(tests, function(x) x[2, "p-value"])
) # all variables are significant additions 
tests


tab_model(lme_complex_ar1, file = 'complexity.HTML')


############# F. PLOTS #############
# Visualise random effects 
#dive duration
(re.effects_dd <- plot_model(lme_dur_ar1, type = "est", 
                         rm.terms = c("substrateCoarse", 'substrateRocky', "sunElevationNight", "bathymetry"),
                         title = '', axis.labels = c('Trawl overlap', 'Gillnet overlap'),
                         vline.color = "#061423", dot.size = 5, line.size = 2, value.size = 5.5, colors = 'black')+
   labs(y = 'Dive duration (s)') +
   theme_bw() +
   theme(text = element_text(size = 16, colour = 'black'), 
         legend.text = element_text(size=14, color = 'black'),
         axis.text = element_text(size = 14, color = 'black'),
         panel.grid.major.y = element_blank()))
#time at bottom
(re.effects_sb <- plot_model(lme_sbtm_ar1, type = "est", 
                          rm.terms = c("substrateCoarse", 'substrateRocky', "sunElevationNight", "bathymetry"),
                          title = '', axis.labels = c('.', '.'),
                          vline.color = "#061423", dot.size = 5, line.size = 2, value.size = 5.5, colors = 'black')+
    labs(y = 'Time at bottom (s)') +
    theme_bw() +
    theme(text = element_text(size = 16, colour = 'black'), 
          legend.text = element_text(size=14, color = 'black'),
          axis.text = element_text(size = 14, color = 'black'),
          axis.text.y = element_text(color = 'white'),
          panel.grid.major.y = element_blank()))
#complexity
(re.effects_c <- plot_model(lme_complex_ar1, type = "est", 
                          rm.terms = c("substrateCoarse", 'substrateRocky', "sunElevationNight", "bathymetry"),
                          title = '', axis.labels = c('Trawl overlap', 'Gillnet overlap'), colors = 'black',
                          vline.color = "#061423", dot.size = 5, line.size = 2)+
    labs(y = 'Complexity (m)') +
    theme_bw() +
    theme(text = element_text(size = 16, colour = 'black'), 
          legend.text = element_text(size=14, color = 'black'),
          axis.text = element_text(size = 14, color = 'black'),
          panel.grid.major.y = element_blank())+
    annotate('text', label = 'p < 0.001', x = 2.16, y = 2.5))
#swim effort 
(re.effects_se <- plot_model(lme_swimeff_ar1, type = "est", 
                          rm.terms = c("substrateCoarse", 'substrateRocky', "sunElevationNight", "bathymetry"),
                          title = "", axis.labels = c('.', '.'),
                          digits = 4, vline.color = "#061423", axis.lim = c(-1.5, 1.5),
                          dot.size = 5, line.size = 2, value.size = 5.5, colors = 'black') + 
    labs(y = 'Swim effort') +
    theme_bw() +
    theme(text = element_text(size = 16, colour = 'black'), 
          legend.text = element_text(size=14, color = 'black'),
          axis.text = element_text(size = 14, color = 'black'),
          axis.text.y = element_text(color = 'white'),
          panel.grid.major.y = element_blank())+
    ylim(-0.01, 0.01))

(fig5 <- ((re.effects_dd + re.effects_sb) / (re.effects_c + re.effects_se)) +
  plot_annotation(tag_levels = 'A'))

#model output tabels 
tab_model(lme_dur_ar1,
          pred.labels =c("(Intercept)", "Gillnet Overlap", "Trawl Overlap", "Coarse Substrate",
                         'Rocky Substrate', 'Night', 'Bathymetry'),
          dv.labels= "Dive Duration", file = 'divedur.HTML')
tab_model(lme_sbtm_ar1,
          pred.labels =c("(Intercept)", "Gillnet Overlap", "Trawl Overlap", "Coarse Substrate",
                         'Rocky Substrate', 'Night', 'Bathymetry'),
          dv.labels= "Time at Bottom", file = 'secsbtm.HTML')
tab_model(lme_complex_ar1,
          pred.labels =c("(Intercept)", "Gillnet Overlap", "Trawl Overlap", "Coarse Substrate",
                         'Rocky Substrate', 'Night', 'Bathymetry'),
          dv.labels= "Complexity", file = 'complexity.HTML')
tab_model(lme_swimeff_ar1,
          pred.labels =c("(Intercept)", "Gillnet Overlap", "Trawl Overlap", "Coarse Substrate",
                         'Rocky Substrate', 'Night', 'Bathymetry'),
          dv.labels= "Swim Effort", file = 'swimeff.HTML')



