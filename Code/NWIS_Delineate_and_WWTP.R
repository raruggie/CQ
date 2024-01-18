# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

library(FactoMineR)
library(factoextra)
library(leaps)
library(glmnet)
library(randomForest)
library(ggcorrplot)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(climateR)
library(ggpubr)
library(ggnewscale)
library(ggforce)
library(readxl)
library(MASS)
library(car)
library(CropScapeR)
library(FedData)
library(streamstats)
library(dataRetrieval)
library(sf)
library(raster)
library(terra)
library(kableExtra)
library(sfheaders)
library(mapview)
library(broom)
library(ggsignif)
library(ggpmisc)
library(segmented)
library(readxl)
library(tmap)
library(readxl)
library(zyp)
library(corrr)
library(tidyverse)

sf_use_s2(TRUE) # for sf

setTimeout(1000) # for streamstats api

meters_to_miles = 1/1609.334

####################### Functions #######################

source("Code/Ryan_functions.R")

####################### Goal of code #######################

# 1) Process the NWIS database for the CQ analysis
# 2) Run through CQ metrics/watershed attribute correlation

####################### Workflow #######################

# read in site metadata and the df.TP_CQ df to get the sites to filter on:

load('Processed_Data/TP.Rdata')

df.NWIS.TP_site_metadata<-read.csv("Raw_Data/df.NWIS.TP_site_metadata.csv", colClasses = c(site_no = "character"))%>%
  filter(site_no %in% df.TP_CQ$site_no)

#### Delinating Watersheds ####

# delinate the sites:

# l.SS_WS.NWIS<-lapply(seq_along(df.NWIS.TP_site_metadata$site_no), \(i) Delineate(df.NWIS.TP_site_metadata$dec_long_va[i], df.NWIS.TP_site_metadata$dec_lat_va[i]))

# names(l.SS_WS.NWIS)<-df.NWIS.TP_site_metadata$site_no

# save(l.SS_WS.NWIS, file = 'C:/PhD/CQ/Downloaded_Data/NWIS.TP.SS_WS.Rdata')

# df.sf.NWIS<-fun.l.SS_WS.to.sfdf(l.SS_WS.NWIS)

# save(df.sf.NWIS, file = 'C:/PhD/CQ/Processed_Data/NWIS.TP.SS_WS_sf.Rdata')

load('C:/PhD/CQ/Processed_Data/NWIS.TP.SS_WS_sf.Rdata')

# which sites didnt work:

# first the ones that didnt make it out of going from SS_WS to df.sf:
  
setdiff(df.NWIS.TP_site_metadata$site_no, df.sf.NWIS$Name)

# then I used the streamstats webpage to get the watersheds for these
# reading them in:

WS_01357500<-st_read('Raw_Data/01357500_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '01357500')
WS_01362200<-st_read('Raw_Data/01362200_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '01362200')
WS_04213500<-st_read('Raw_Data/04213500_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '04213500')

# now combine these and add to df.sf.NWIS:

temp<-bind_rows(WS_01357500, WS_01362200, WS_04213500)

mapview(temp)

df.sf.NWIS<-bind_rows(df.sf.NWIS, temp)

# calculate DA using vect:

vect.NWIS<-vect(df.sf.NWIS)

vect.NWIS$area_KM2<-expanse(vect.NWIS, unit="km")

# add DA in mi2 to df.sf:

df.sf.NWIS$area_sqmi<-expanse(vect.NWIS, unit="km")*0.386102

# add NWIS tabulated area to df sf: first download the metadata for the sites using readNWISsite, then merge drain_area_va column to sfdf:

temp<-readNWISsite(df.sf.NWIS$Name)

df.sf.NWIS<-left_join(df.sf.NWIS, temp[,c(2,30)], by = c('Name'='site_no'))

# calcuate the percent error of the delination:

df.sf.NWIS$Delination_Error<-(df.sf.NWIS$area_sqmi-df.sf.NWIS$drain_area_va)/df.sf.NWIS$drain_area_va

# I want to look at the delimaitons that are on the cusp of not working:

temp<-df.sf.NWIS%>%
  filter(abs(Delination_Error)>.02)

mapview(temp)

temp$Name

# then there are these, the oneswith high delineation error.
# going back to streamstats to get these:

WS_01362500<-st_read('Raw_Data/01362500_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '01362500')
WS_04231600<-st_read('Raw_Data/04231600_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '04231600')
WS_04249000<-st_read('Raw_Data/04249000_SS_online_WS_shapefile/globalwatershed.shp')%>%mutate(Name = '04249000')

# now combine these and add to df.sf.NWIS. to do this:

# reload in df.sf.NWIS:

load('C:/PhD/CQ/Processed_Data/NWIS.TP.SS_WS_sf.Rdata')

# remove these three:

df.sf.NWIS<-filter(df.sf.NWIS, !Name %in% temp$Name)

# add these three plus the three from before into temp:

temp<-bind_rows(WS_01362500, WS_04231600, WS_04249000, WS_01357500, WS_01362200, WS_04213500)

# now run through the workflow above again:

df.sf.NWIS<-bind_rows(df.sf.NWIS, temp)
vect.NWIS<-vect(df.sf.NWIS)
vect.NWIS$area_KM2<-expanse(vect.NWIS, unit="km")
df.sf.NWIS$area_sqmi<-expanse(vect.NWIS, unit="km")*0.386102
temp<-readNWISsite(df.sf.NWIS$Name)
df.sf.NWIS<-left_join(df.sf.NWIS, temp[,c(2,30)], by = c('Name'='site_no'))
df.sf.NWIS$Delination_Error<-(df.sf.NWIS$area_sqmi-df.sf.NWIS$drain_area_va)/df.sf.NWIS$drain_area_va

# now check the delineation errors:

temp<-df.sf.NWIS%>%
  filter(abs(Delination_Error)>.02)

# good, there are none!

# export watershed shapefile sf df:

save(df.sf.NWIS, file = 'Processed_Data/NWIS_Watershed_Shapefiles.Rdata')

#### finding WWTP in watersheds ####

# read in shapefile of NYS WWTP:

wwtp<-st_read('Raw_Data/NYS_Wastewater_Facility_shapefile/Wastewater_Facility.shp')

# create map of watersheds and WWTP:

mapview(df.sf.NWIS)+mapview(wwtp)

# find intersection of WWTP and watersheds:

x<-st_intersection(wwtp, df.sf.NWIS)

mapview(df.sf.NWIS)+mapview(x)

# looks good
# write out to csv:

write.csv(x[,c(1:6)], row.names = FALSE, file = 'Processed_Data/WWTP_in_NWIS_WS.csv')

