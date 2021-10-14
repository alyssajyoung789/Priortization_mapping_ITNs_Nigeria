# --------------------------------------------------------------------------
# CREATING MAPS FOR NET TARGET PUBLICATION --------------------------------
# Author: Alyssa Young, Tulane SPHTM --------------------------------------
# Date last modified: 10/14/2021 -------------------------------------------
# --------------------------------------------------------------------------


# Load packages
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
library (tmap)

# Explore color palette ---------------------------------
palette_explorer()
magma <- viridisLite::magma(6)
plasma <- viridisLite::plasma(6)
inferno <- viridisLite::inferno(6)
viridis <- viridisLite:: viridis(6)
cividis<-viridisLite:: cividis(6)

# Get admin boundaries ----------------------------------------------------

adm0.nga <- st_read(dsn = "C:\\Users\\Alyssa\\OneDrive\\Desktop\\Malaria Consortium\\Nigeria data\\Geospatial\\gadm36_NGA_0.shp")
adm1.nga <- st_read(dsn = "C:\\Users\\Alyssa\\OneDrive\\Desktop\\Malaria Consortium\\Nigeria data\\Geospatial\\gadm36_NGA_1.shp")
adm2.nga <- st_read(dsn = "C:\\Users\\Alyssa\\OneDrive\\Desktop\\Malaria Consortium\\Nigeria data\\Geospatial\\gadm36_NGA_2.shp")

plot(adm0.nga)
plot(adm1.nga)
plot(adm2.nga)
tm_shape(adm2.nga) + tm_polygons("NAME_2") + tm_legend(show=FALSE)

# View adm2.nga as data frame
adm2.nga.df <- as.data.frame(adm2.nga)
adm1.nga.df <- as.data.frame(adm1.nga)


# View adm2.nga.df
View(adm2.nga.df)

### reading in dataset to map
NGA_PUB <-read.csv("C:\\Users\\Alyssa\\OneDrive\\Desktop\\Malaria Consortium\\Nigeria data\\NGA_maps_pub.csv")
NGA_PUB.df <- as.data.frame(NGA_PUB)
NGA_PUB.spdf <- merge (adm2.nga, NGA_PUB.df)


# --------------------------------------------------------------------------
# OBJECTIVE 1: create one figure of all input map layers used to -------
#create final prioritization map-------------------------------------------
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# 1) Create mean annual rainfall map ---------------------------------------
# --------------------------------------------------------------------------


RF_map <- tm_shape(NGA_PUB.spdf) + 
  tm_fill("rf_class", title = "\nMean \nannual \nrainfall (mm)",
          n = 4, style = "cat", 
          palette = "Greens", labels = c("\U2264 91", " 92 - 147", "148 - 207", "> 207"), 
          lwd =0.2, border.col = "grey40") + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.2, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "a", title.size= 1, title.snap.to.legend= TRUE, frame = FALSE, legend.outside = TRUE, 
            legend.position = c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)


# --------------------------------------------------------------------------
# 2) Create temperature suitability map ------------------------------------
# --------------------------------------------------------------------------
# create custom color palette
YlOrBr <- get_brewer_pal("YlOrBr", n = 4)


temp_suit_map <- tm_shape(NGA_PUB.spdf) +
  tm_fill("temp_suit_class", title = "\nP.f.\ntemperature\nsuitability\nindex",
          n = 4, style = "cat",
          palette = YlOrBr, labels = c("> 0.680", "0.571 - 0.680", "0.454 - 0.570", "\U2264  0.453"),
          lwd =0.2, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.2, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "b", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)


# --------------------------------------------------------------------------
# 3) Create No. years since last ITN distribution  map ---------------------
# --------------------------------------------------------------------------

ITN_years_map <- tm_shape(NGA_PUB.spdf) + 
  tm_fill("ITN_dist_class", title = "\nYears since\n last ITN\n distribution",
          n = 4, style = "cat", 
          palette = "Blues", labels = c("< 3 years", "3 years",
                                        "5 years", "\u2265 6 years"), 
          lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "c", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)

# --------------------------------------------------------------------------
# 4) Create PBO map --------------------------------------------------------
# --------------------------------------------------------------------------

BuGn <- get_brewer_pal("GnBu", n = 2, contrast = c(0, 0.42))


pbo_map<- tm_shape(NGA_PUB.spdf) + 
  tm_fill("PBO_class", title = "\nPBO nets \ndistributed\nin 2019  ", n=2, style = "cat", 
          palette = "GnBu", labels = c("Yes", "No"), 
          lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "d", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)

# --------------------------------------------------------------------------
# 5a) Create ITN access map v1-----------------------------------------------
# --------------------------------------------------------------------------


RdPu <- get_brewer_pal("RdPu", n = 4, contrast = c(0, 0.42))

ITN_INLA_map <- tm_shape(NGA_PUB.spdf) + 
  tm_fill("ITN_access_INLA_mclass", title = "\nITN Access v1\n(INLA)", n=4, style = "cat", 
          palette = "RdPu", labels = c("\u2265 49%", "35 - 48%", "23 - 34%", "<23%"),
          lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "j", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)

# --------------------------------------------------------------------------
# 5b) Create ITN availbility map v2-----------------------------------------
# --------------------------------------------------------------------------


RdPu <- get_brewer_pal("RdPu", n = 4, contrast = c(0, 0.42))

ITN_DHS_map <- tm_shape(NGA_PUB.spdf) + 
  tm_fill("ITN_access_DHS_18_state_precent_class", title = "\nITN Access v2\n(DHS)", n=4, style = "cat", 
          palette = "RdPu", labels = c("> 45.2 %", "33 - 45.2%", "19 - 32.9%", "<18%"),
          lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") + 
  tm_scale_bar(position = c(0.55, 0.01)) + tm_compass(position= c(0.8, 0.2)) +
  tm_layout(title= "k", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)





# --------------------------------------------------------------------------
# 6) Create SMC map --------------------------------------------------------
# --------------------------------------------------------------------------
ltred <- get_brewer_pal("-Reds", n = 3, contrast = c(0.04, 0.44))

SMC_map <-tm_shape(NGA_PUB.spdf) + 
  tm_fill("smc_class", title = "SMC LGAs", n=2, style = "cat", 
          palette = ltred, labels = c("Yes", "No"), 
          lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "e", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)

# --------------------------------------------------------------------------
# 7) Create Built up presence map (SMOD) ----------------------------------
# --------------------------------------------------------------------------

# create custom color palette
YlGnBu <- get_brewer_pal("YlGnBu", n = 4, contrast = c(0, 0.42))

SMOD_map <-tm_shape(NGA_PUB.spdf) +
  tm_fill("Built_class", title = "\nBuilt up\npresence\nindex", n = 4, style = "cat",
          palette = YlGnBu, labels = c(" \U2264 0.0085", "0.051 -  0.0084", "0.76 - 0.05", "> 0.76"),
  lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "f", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)


# --------------------------------------------------------------------------
# 8) Create IDP map --------------------------------------------------------
# --------------------------------------------------------------------------
IDP_map <- tm_shape(NGA_PUB.spdf) + 
  tm_fill("IDP_cat_binary", title = "Internally\ndisplaced\npopulations \n  ", n=2, style = "cat", 
          palette = "Purples", labels = c("No", "Yes"), 
          lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "g", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)

# --------------------------------------------------------------------------
# 9) Create Malaria prevalence map -----------------------------------------
# --------------------------------------------------------------------------
Prev_map <-tm_shape(NGA_PUB.spdf) +
  tm_fill("Prev_DHS_18_state_class", title = "Malaria \nPrevalence \n6-59 mo", n = 5, style = "cat",
          palette = "Reds", labels = c(" <12", "13 - 22%", "23 - 32%", "33 - 42%", ">42%"),
          lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) + 
  tm_layout(title= "h", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)

# --------------------------------------------------------------------------
# 10) Create SES map --------------------------------------------------------
# -------------------------------------------------------------------------- 

BuGn <- get_brewer_pal("GnBu", n = 4, contrast = c(0, 0.50))

SES_map <-tm_shape(NGA_PUB.spdf) +
  tm_fill("Perc_pop_lowest_wealth_quintile_DHS_2018_class", title = "% Population\nlowest\nquintile", n = 4, style = "cat",
          palette = BuGn, labels = c(" <8% ", "8 - 24.5%", "24.6 - 40.8%", "> 40.8%"),
          lwd =0.5, border.col = "grey40")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 0.5, col = "grey40") +
  tm_scale_bar(position = c(0.55, 0.01)) +  
  tm_layout(title= "i", title.size= 1, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "gray39", lwd = 1)


# --------------------------------------------------------------------------
# OBJECTIVE 2: create one figure of both priortization scheme maps
#  The first map features spatially interpolated ITN access layer as an input,
# the other uses DHS state-reported ITN access layer as an input
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# 12a) Create Weighted Priortization/Malaria Risk map with INLA access layer
# --------------------------------------------------------------------------

Final_prioritization_map_INLA <-tm_shape(NGA_PUB.spdf) +
  tm_fill("Final_Risk_score_INLA_access_class", title = "Nigeria \nLGA Priority \nRankings (v1)", n = 6, style = "cat",
          palette = "Blues", labels = c("Low", "Moderate Low", "Moderate", "Moderate High", "High", "Extremely High"), 
          lwd =0.5, border.col = "gray39")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 1, col = "gray39") +
  tm_scale_bar(position = c(0.55, 0.01)) + tm_compass(position= c(0.8, 0.2)) +
  tm_layout(title= "a", title.size= 1.5, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "grey14", lwd = 1)


# --------------------------------------------------------------------------
# 12b) Create Weighted Priortization/Malaria Risk map with DHS access layer-
# --------------------------------------------------------------------------

# with DHS jenks
Final_prioritization_map_DHS <-tm_shape(NGA_PUB.spdf) +
  tm_fill("Final_Risk_Score_DHS_access_class", title = "Nigeria \nLGA Priority \nRankings (v2)", n = 6, style = "cat",
          palette = "Blues", labels = c("Low", "Moderate Low", "Moderate", "Moderate High", "High", "Extremely High"), 
          lwd =0.5, border.col = "gray39")  + 
  tm_shape(adm2.nga)+ tm_borders (lwd = 1, col = "gray39") +
  tm_scale_bar(position = c(0.55, 0.01)) + tm_compass(position= c(0.8, 0.2)) +
  tm_layout(title= "b", title.size= 1.5, title.snap.to.legend= TRUE,frame = FALSE, legend.outside = TRUE,
            legend.position =c("RIGHT","BOTTOM"), legend.title.size = 1, legend.text.size= 0.65, legend.hist.size = 0.65) +
  tm_shape(adm1.nga) + tm_borders(col = "grey14", lwd = 1)

# --------------------------------------------------------------------------
# Make tmap_arrange figure -------------------------------------------------
# --------------------------------------------------------------------------

tmap_arrange_1 <- tmap_arrange(RF_map, 
                               temp_suit_map, 
                               ITN_years_map,
                               pbo_map,
                               SMC_map,
                               SMOD_map,
                               IDP_map,
                               Prev_map,
                               SES_map,
                               ITN_INLA_map,
                               Final_prioritization_map_DHS,
                               Final_prioritization_map_INLA,
                               ncol = 2,
                               outer.margins = 0)

# Save as tmap arrange (multiple maps in same figure)

tmap_save(tmap_arrange_1, 
          "final_series.png", width =7,
          height=7, units ='in', asp = 0)





