# ------------------------------------------------------------------------------------------
# CREATING ITN DISTRIBUTION PRIORTIZATION SCHEME FOR NIGERIA
# Objective 1 :Create classes for each input factor using Jenks natural breaks or presence/absence of intervention
# Objective 2: Assign weights to input factors using AHP methodology 
# Author(s): Alyssa Young, Will Eaton, Tulane SPHTM 
# Date last modified: 10/14/2021 
# ------------------------------------------------------------------------------------------

# load libraries
library(ahpsurvey)
library(raster)
library(sp)
library(tmap)
library(tmaptools)
library(BAMMtools)
library(readxl)
library(ggplot2)
library(malariaAtlas)
library(exactextractr)
library(maptools)
library(raster)
library(sf)
library("spDataLarge")
library(dplyr)
library(GADMTools)

#library(exactextractr)

library(maptools)
library(readr)
library(foreign)   ; library(tsModel) ; library("lmtest") ; library("Epi")
library("splines") ; library("vcd")
library(reshape2)  ; library(hablar)
library(tidyr)     
library (viridis)
library(data.table)
library(forecast)  ; library(MASS)
library(tseries)   ; library(scales)
library(tsModel)   ; library(extrafont)
library(lmtest)    ; library(tidyverse)
library(stargazer) ; library(RColorBrewer)
library(readxl)    ; library(olsrr)
library(Hmisc)
library(MASS)
library(dplyr)
library(devEMF)
library(padr)
library(zoo)
library(tidyverse)
library(naniar)
library(GGally)
library(cartogram)
library(mgcv)
library(BAMMtools)
library(rgdal)
library(leaflet)
library(arsenal)
require(data.table)
require(RCurl)
require(R.utils)
require(gdalUtils)
require(parallel)

## load dataset; data organized by admin 2 name factor weights, factor scores, and final priortizaiton score per LGA 
## as variables 
NGA_PUB <-read.csv("C:\\Users\\Alyssa\\OneDrive\\Desktop\\Malaria Consortium\\NGA_pub_maps_AHP_wt_10_21.csv")

#-----------------------------------------------------------------------------------------------
# Create Jenks classes based on distribution of data for each input factor/variable-------------
# ----------------------------------------------------------------------------------------------
# With natural breaks classification (Jenks) Natural Breaks Jenks, 
# classes are based on natural groupings inherent in the data. 
# Class breaks are identified that best group similar values and 
# that maximize the differences between classes. The features are 
# divided into classes whose boundaries are set where there are 
# relatively big differences in the data values.
# 
# Natural breaks are data-specific classifications and not useful for 
# comparing multiple maps built from different underlying information.

# 1.) State-level ITN access (ITN Access) ------------------------------------------------------

getJenksBreaks(NGA_PUB$ITN_access_INLA, 5, subset = NULL)  
# [1] 0.009014783 0.233290972 0.354107688 0.493038245 0.792025051

NGA_PUB$ITN_access_INLA_class[NGA_PUB$ITN_access_INLA < 23.0] <- 4
NGA_PUB$ITN_access_INLA_class[NGA_PUB$ITN_access_INLA <35.0 & NGA_PUB$ITN_access_INLA >= 23.0 ]  <- 3
NGA_PUB$ITN_access_INLA_class[NGA_PUB$ITN_access_INLA < 48.0 & NGA_PUB$ITN_access_INLA >= 35.0] <- 2
NGA_PUB$ITN_access_INLA_class[NGA_PUB$ITN_access_INLA >= 48.0] <- 1

# 2.) Mean annual rainfall (RFE) ----------------------------------------------------------------------------
getJenksBreaks(NGA_PUB$mean_an_rf, 5, subset = NULL)  
#[1] -53.22974  92.29017 146.79149 206.53586 300.44352

NGA_PUB$rf_class[NGA_PUB$mean_an_rf > 207.0] <- 4
NGA_PUB$rf_class[NGA_PUB$mean_an_rf <= 207 & NGA_PUB$ITN_access_INLA > 147.0]  <- 3
NGA_PUB$rf_class[NGA_PUB$mean_an_rf <= 147.0 & NGA_PUB$ITN_access_INLA >= 92.0] <- 2
NGA_PUB$rf_class[NGA_PUB$mean_an_rf <= 91.0] <- 1

# 3.) Built up area presence (SMOD): proxy for rural/urban designation ------------------------------------------------------------
getJenksBreaks(NGA_PUB$SMOD, 5, subset = NULL)

NGA_PUB$Built_class[NGA_PUB$SMOD <= 0.0084] <- 4
NGA_PUB$Built_class[NGA_PUB$SMOD <= 0.051 & NGA_PUB$SMOD >= 0.0085]  <- 3
NGA_PUB$Built_class[NGA_PUB$SMOD <= 0.76 & NGA_PUB$SMOD >= 0.05] <- 2
NGA_PUB$Built_class[NGA_PUB$SMOD > 0.76] <- 1

# 4.) Plasmodium falciparum temperature suitability index (TSI) --------------------------------------------------------------

# Show malaria atlas project data available
listData('raster', printed = TRUE)
listRaster(printed = TRUE)

# Download malaria atlas project temp suitability raster layer 
NGA_shp <- getShp(ISO = "NGA", admin_level = "admin0")

# (From https://rdrr.io/cran/malariaAtlas/man/getRaster.html)
nga_temp_suit_raster <- getRaster(surface = "Plasmodium falciparum Temperature Suitability", shp = NGA_shp)
# Visualize raster in plot
autoplot_MAPraster(nga_temp_suit_raster)

getJenksBreaks(NGA_PUB$mean_temp_suit, 5, subset = NULL)
# [1] 0.3015740 0.4532296 0.5710385 0.6795321 0.7636224

NGA_PUB$temp_suit_class[NGA_PUB$mean_temp_suit > 0.680] <- 4
NGA_PUB$temp_suit_class[NGA_PUB$mean_temp_suit <= 0.680 & NGA_PUB$SMOD > 0.570]  <- 3
NGA_PUB$temp_suit_class[NGA_PUB$mean_temp_suit <= 0.570 & NGA_PUB$SMOD >= 0.454] <- 2
NGA_PUB$temp_suit_class[NGA_PUB$mean_temp_suit < 0.454] <- 1


# 5.) PBO nets previously distributed (PBO)
# 6.) Number of years since last ITN distribtuion (NoYrsITN)
# 7.) Presence of internally displaced persons (IDPs)
# 8.) DHS-reported state level microscopy malaria prevalence among children under 5 (Prev)


# 9.) Socioeconomic status/Mean Wealth Index (MWI): Percent pop in lowest quintile per state ---------------
getJenksBreaks(NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018,5, subset = NULL)  
# [1]  0.0  7.9 24.5 40.8 63.2
#Assign Percent pop in lowest quintile per state
NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018_class[NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018 > 40.8] <- 4     # classify as rank 4
NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018_class[NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018 <= 40.8 & NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018 > 24.5 ] <- 3   # classify as rank 3
NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018_class[NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018 <= 24.5 & NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018 > 7.9 ] <- 2    # classify as rank 2
NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018_class[NGA_PUB$Perc_pop_lowest_wealth_quintile_DHS_2018 <= 7.9] <- 1 # classify as rank 1


#-----------------------------------------------------------------------------------------------
# Obtain weights, using AHP (canned approach) for each input factor/variable--------------------
# ----------------------------------------------------------------------------------------------
# Code source: https://cran.r-project.org/web/packages/ahpsurvey/vignettes/my-vignette.html

AHP.NGA <-read.csv("C:\\Users\\Alyssa\\OneDrive\\Desktop\\Malaria Consortium\\AHP_translated_survey_response.csv")
atts<- c("RFE", "PfTSI", "prev", "ITNAccess", "NoyrsITN", "PBO", "UrbRur", "SMC", "MWI", "IDPs")

canned <- ahp(df = AHP.NGA, 
              atts = c('RFE', 'PfTSI', 'prev', 'ITNAccess', 'NoyrsITN', 'PBO', 'UrbRur', 'SMC', 'MWI', 'IDPs'), 
              negconvert = FALSE, 
              reciprocal = TRUE,
              method = 'arithmetic', 
              aggmethod = "arithmetic", 
              qt = 0.2,
              censorcr = 0.1, # setting up code so that any observations that would make the CR value over than 0.1 are dropped
              agg = TRUE)

head(canned$indpref)
canned$aggpref

# no observations censored

#            AggPref  SD.AggPref
# RFE       0.16185757 0.09481153
# PfTSI     0.08861157 0.07166681
# prev      0.06793341 0.04615549
# ITNAccess 0.14756037 0.06526230
# NoyrsITN  0.18081495 0.08176559
# PBO       0.09556061 0.07052058
# UrbRur    0.08734026 0.07117833
# SMC       0.08138323 0.05804888
# MWI       0.06174643 0.02792138
# IDPs      0.02719159 0.00903272


#### CREATE CLASSES BASED ON FINAL PRIORITIZATION SCORES USING BOTH ITN ACCESS LAYERS 

# load data
NGA_PUB <-read.csv("C:\\Users\\Alyssa\\OneDrive\\Desktop\\Malaria Consortium\\NGA_pub_maps_AHP_wt_10_21.csv")

# 1a. )Final risk categorization version 1 (with INLA ITN Access scores) -----------------------------------------------------------

getJenksBreaks(NGA_PUB$Final_Risk_score_INLA_access, 5, subset = NULL) 
# [1] 1.37 1.82 2.16 2.54 3.06

NGA_PUB$Final_Risk_score_INLA_access_class [NGA_PUB$Final_Risk_score_INLA_access >= 3.06] <- 6 # Extremely High
NGA_PUB$Final_Risk_score_INLA_access_class [NGA_PUB$Final_Risk_score_INLA_access < 3.06 & NGA_PUB$Final_Risk_score_INLA_access > 2.54] <- 5 # High
NGA_PUB$Final_Risk_score_INLA_access_class [NGA_PUB$Final_Risk_score_INLA_access <= 2.54 & NGA_PUB$Final_Risk_score_INLA_access > 2.16] <- 4 # Moderate High
NGA_PUB$Final_Risk_score_INLA_access_class [NGA_PUB$Final_Risk_score_INLA_access <= 2.16 & NGA_PUB$Final_Risk_score_INLA_access > 1.82] <- 3 # Moderate
NGA_PUB$Final_Risk_score_INLA_access_class [NGA_PUB$Final_Risk_score_INLA_access <= 1.82 & NGA_PUB$Final_Risk_score_INLA_access > 1.37] <- 2 # Moderate Low
NGA_PUB$Final_Risk_score_INLA_access_class [NGA_PUB$Final_Risk_score_INLA_access <= 1.37 ] <- 1 # Low


# 1b.) Final risk categorization version 2.1 ( with Sate-level DHS ITN Access scores)----------------------------------------------------------
getJenksBreaks(NGA_PUB$Final_Risk_Score_DHS_access, 5, subset = NULL) 
#  [1] 1.19 1.81 2.16 2.49 2.91

NGA_PUB$Final_Risk_Score_DHS_access_class [NGA_PUB$Final_Risk_Score_DHS_access >= 2.91] <- 6 # Extremely High
NGA_PUB$Final_Risk_Score_DHS_access_class [NGA_PUB$Final_Risk_Score_DHS_access < 2.91 & NGA_PUB$Final_Risk_Score_DHS_access > 2.49] <- 5 # High
NGA_PUB$Final_Risk_Score_DHS_access_class [NGA_PUB$Final_Risk_Score_DHS_access <= 2.49 & NGA_PUB$Final_Risk_Score_DHS_access > 2.16] <- 4 # Moderate High
NGA_PUB$Final_Risk_Score_DHS_access_class [NGA_PUB$Final_Risk_Score_DHS_access <= 2.16 & NGA_PUB$Final_Risk_Score_DHS_access > 1.81] <- 3 # Moderate
NGA_PUB$Final_Risk_Score_DHS_access_class [NGA_PUB$Final_Risk_Score_DHS_access <= 1.81 & NGA_PUB$Final_Risk_Score_DHS_access > 1.19] <- 2 # Moderate Low
NGA_PUB$Final_Risk_Score_DHS_access_class [NGA_PUB$Final_Risk_Score_DHS_access <= 1.19] <- 1 # Low
