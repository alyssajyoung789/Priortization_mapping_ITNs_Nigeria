# ------------------------------------------------------------------------------------------
# CREATING ITN DISTRIBUTION PRIORTIZATION SCHEME FOR NIGERIA
# Objective 1 :Create classes for each input factor using Jenks natural breaks or presence/absence of intervention
# Objective 2: Assign weights to input factors using AHP methodology 
# Author: Alyssa Young, Tulane SPHTM 
# Date last modified: 10/14/2021 
# ------------------------------------------------------------------------------------------

# load libraries
library(ahpsurvey)
library(raster)
library(sp)
library(BAMMtools)
library(exactextractr)

## load dataset; data organized by admin 2 name factor weights, factor scores, and final priortizaiton score per LGA 
## as variables 
NGA_PUB <-read.csv("C:\\Users\\Alyssa\\OneDrive\\Desktop\\Malaria Consortium\\Nigeria data\\NGA_maps_pub.csv")


#-----------------------------------------------------------------------------------------------
# Create Jenks classes based on distribution of data for each input factor/variable-------------
# ----------------------------------------------------------------------------------------------

# 1.) State-level ITN access (ITN Access) ------------------------------------------------------
# 2.) Mean annual rainfall (RFE)
# 3.) Built up area presence (SMOD): proxy for rural/urban designation
# 4.) Plasmodium falciparum temperature suitability indiex (TSI)
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
