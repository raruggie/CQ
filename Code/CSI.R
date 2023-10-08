# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

library(climateR)
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

####################### Functions #######################

source("C:/PhD/CQ/Code/Ryan_functions.R")

####################### Goal of code #######################

# Process the CSI database for the CQ analysis

####################### Workflow #######################

# Like DOW, flow data will need to be paired up for the CSI sampling sites
# however, there are only two watersheds. 

# I'm getting fustrated thinking what Chuck and Steve are going to say when we look at the NWIS and 
# then the CSI CQ curves. Its like comparing apples and oranges. We have a lot of NWIS sites, so why not just use those in the analysis since we trust the flow data?

# Also, the CSI sites have so many fewer samples, is it bias to compare them with the NWIS sites?

# I want to at least start the CQ analysis for Chuck and Steve
# I can always come back to this code and add more sites if we want.

# Below is theworkflow when the threshold was n > 50 samples: (I just liked 99 for DOW and 50 for CSI, no reasoningotherthan how Iwas feeling when I startedworking onthe datasets)

# the code is adapted from C:/PhD/Research/Code/Data_Inventory/Query_and_export_CSI_data.R

# Read in Raw CSI data downloaded from their website and saved as csv: to do this:

setwd("C:/PhD/Research/Data/Data_inventory/CSI")

# import the raw C data as a csv:

df.CSI<-read.csv("CSI-all_monitoring_regions.csv")

# format the dataframe: to do this:

# remove whitespace from monitoring set and location columns:

df.CSI$Monitoring.Set<-str_trim(df.CSI$Monitoring.Set, side="both")

df.CSI$Monitoring.Location<-str_trim(df.CSI$Monitoring.Location, side = "both")

# convert the date: 

df.CSI$Date<-as.Date(df.CSI$Date, "%b %d, %Y")

# create a unique site ID column that combines set and location name:

df.CSI$Name<-paste(df.CSI$Monitoring.Set, "-", df.CSI$Monitoring.Location)

# pair Raw CSI WQ data to the sites Lat/Long Data. Each sub watersheds sites' lat/longs were saved in a nested file strucuter, 
# To do this: 

# extract the names of the csv files in the nested folder structure into a list of 2 characters objects, on of the watershed and one of the monitoring set (which is the name of CSV in the monitring set folder). The observations in the csv moniotring set are the lat longs of the monitiring loicaitons: :

filenames <- list.files("C:/PhD/Research/Data/Data_inventory/CSI/lat_long", recursive=TRUE)%>%
  str_split(., "/")

# combine the monitoring set and monitoring location (i.e. file name) into a 2 columndataframe:

filenames<-as.data.frame(do.call(rbind, filenames))

# set the column names of the dataframe:

names(filenames)<-c("Watershed", "Monitoring_Set")

# get the unique watersheds in a vector:

watersheds<-unique(filenames$Watershed)

# create a vector of the unique monitoring sets by removing the .csv ending in the file name:

Monitoring_set_sites<-sort(substr(unique(filenames$Monitoring_Set), 1,nchar(unique(filenames$Monitoring_Set))-4))

# create an empty dataframe of the site (set+location) and its lat/long

cn<-c("Name", "Latitude", "Longitude")
df.CSI_sites<-setNames(data.frame(matrix(ncol = length(cn))), cn)

# test iteration for loop below:

ws<-watersheds[5]

# loop through the watersheds (folders) and then the monitoring sets (csv files) to extract the location name and the lat longs and add to dataframe:

for (ws in watersheds){
  
  # set working director to the watershed folder. Will then loop through this folders files in the nested for loop below.
  
  setwd(paste0("C:/PhD/Research/Data/Data_inventory/CSI/lat_long/", ws))
  
  # create a vector of the csv file names in the folder (i.e. a vector of the monitoring sets):
  set <- list.files(pattern="*.csv")
  
  # remove any whitespace on the file name as well as the .csv ending:
  
  set_names <- str_trim(substr(set, 1, nchar(set)-4), side="both")
  
  # test interation for loop below
  
  # i<-3
  
  # loop through the monitoring sets (files in the watershed folder):
  
  for (i in 1:length(set)){
    
    # monitoring location is read in as a csv using the file name in the set vector created in the first loop above
    # and only the location name,lat and long columns are selected:
    
    location<-read.csv(set[i])%>%
      select(Name, Latitude, Longitude)
    
    # remove leading and trailing whitespace:
    
    location$Name<-str_trim(location$Name, side="both")
    
    # drop locations that dont have observations for latitude (nothing we can do about these):
    
    location<-location[!is.na(location$Latitude),]
    
    # set the name of the locaiton to the monitoring set + monitoring location:
    
    location$Name<-paste(set_names[i], '-', location$Name)
    
    # add this locaiton and its lat/long as an observation in empty the dataframe created before the loop:
    
    df.CSI_sites<-rbind(df.CSI_sites, location)
  }
}

# then add the lat longs to the raw C dataframe:

df.CSI<-left_join(df.CSI, df.CSI_sites, by = 'Name')

# now filter Sites for TP (also reformat the date column):

df.CSI_TP<-df.CSI%>%
  mutate(Date = as.Date(Date, "%Y-%m-%d"))%>%
  filter(grepl('phos', Analyte, ignore.case = T))%>%
  filter(grepl('total', Analyte, ignore.case = T))%>%
  filter(!grepl('dis', Analyte, ignore.case = T))

# note that these duplicate entires may have meaningful data (see next step).
# note that there were a few sites that I could not get the lat longs for (see above loop)
# note that the lat longs were added to the the CSI WQ data download using another code called C:\PhD\CSI\Code\CSI_sites.R

# find unique sites that have lat and long data and have over 50 samples 

df.CSI_TP_sites<-df.CSI_TP%>%
  filter_at(vars(Latitude, Longitude),all_vars(!is.na(.)))%>%
  mutate(ID = paste(Monitoring.Set, Monitoring.Location), .after = 1)%>%
  group_by(ID)%>%
  mutate(n_samples = n(), min_date = min(Date), max_date = max(Date), months = elapsed_months(max(Date),min(Date)))%>%
  distinct(ID, Latitude, Longitude, .keep_all = T)%>%
  arrange(desc(n_samples))%>%
  filter(n_samples >20)

# look at the map of sites:
# first save the sample location points as sf df:

map.CSI.TP_sites<-df.CSI_TP_sites%>%
  as.data.frame(row.names = 1:nrow(.))%>%
  st_as_sf(.,coords=c('Longitude', 'Latitude'), crs = 4326)

mapview(map.CSI.TP_sites, zcol = 'Monitoring.Set')

# looking at this map, there are a limited number of monitoring sets.
# I like the grouping idea of the set so I will assign a surrogate USGS gauge to the set.
# then DA scaling will be used to get flows at each monitoring location in that set.

# Note that the surrogate gauge may not have enough flow data to cover the time period of the sampling
# of each of the monitoring locations. So, OLS or DA scaling maybe be needed to extend the record of the surrogate gauges

# To do this:

# find the start date for flow that I need

sDt<-as.character(min(df.CSI_TP_sites$min_date))

# download NYS USGS flow gauges with data going back this far, and reformat the resulting dataframe of just site and lat long for use in df_dist workflow:

df.CSI.Surro_gauges<-whatNWISsites(statecode = 'NY', parameterCd="00060", startDt=sDt)%>%
  mutate(Site = paste(site_no, station_nm), .before = 1)%>%
  dplyr::select(c(1,3,4,6,7))%>%
  rename(Latitude = dec_lat_va, Longitude = dec_long_va)

# find the surrogate gauge for the monitioring set. To do this:

# the surrogate gauges for the Monitoring Sets on Cayuga Lake were determined by looking at the most downstream sampling site, which had the best chance of being near USGS guage:

# In Cayuga, to do this:

# the following locations are the best for being furthest downstream on their respective tribs:

keep<-c('Virgil Creek Virgil Creek Freeville', 	'Six Mile Creek Plain Street', 'Cayuga Inlet Route 89', "Fall Creek Cayuga Street Bridge", 'Salmon Creek Salmon Creek Mouth', 'Taughannock Creek Taughannock Falls','Trumansburg Creek Camp Barton')

# setting up a dataframe for distance workflow. To do this:filter the CSI TP site to those in keep and just select the name and lat long columns:

CSI_keep<-df.CSI_TP_sites[df.CSI_TP_sites$ID %in% keep,c(2,12,13)]

# rename the site column to x.Site so that Site in the df.CSI.Surro_gauges df can stand alone in the below pipe:

CSI_keep<-CSI_keep%>%rename(x.Site = ID)

# determine all possible combinaitons of the surrogate gauges (quiered above) and CSI locations (the 7in keep). To do this:

# expand() is used here to create all posible combinaitons of the USGS gauges and the 7 CSI downstream sites:
# then left join the lat longs of the CSI sites (already have the gauge lat longs)

df.dist<- df.CSI.Surro_gauges%>%
  group_by(Site, Latitude, Longitude)%>%
  expand(x.Site = CSI_keep$x.Site)%>%
  left_join(.,CSI_keep, by = 'x.Site')%>%
  dplyr::select(c(4,5,6,1,2,3))%>%
  arrange(x.Site)

# rename the columns:

names(df.dist)<-c("x.Site","Latitude.x","Longitude.x","y.Site","Latitude.y","Longitude.y")

# add a distance column (distance between every possible combinaiton of USGS gauge and CSI site) by calling distHaversine (vectorized) on each pair:

df.dist$dist_meters <- geosphere::distHaversine(df.dist[3:2], df.dist[6:5]) # units of meters

# find closest using filter on the distance column:

df.dist<-df.dist%>%
  group_by(x.Site)%>%
  mutate(closest = min(dist_meters))%>%
  filter(closest == dist_meters)

# there are 7 CSI sites and 4 gauges

# add the USGS gauge site_no to the df.dist:

df.dist<-left_join(df.dist, df.CSI.Surro_gauges[,1:2], by = c('y.Site'='Site'))

# get the unique sites from this df:

x<-unique(df.dist$site_no)

# download the raw daily flow data for these sites:

df.CSI.Surro_gauges_Q<-readNWISdv(x, '00060')

# pair these flows up with the sample data:
# first need to merge the monitoring and set back into df.dist:

df.dist<-left_join(df.dist, df.CSI_TP_sites[,2:4], by = c('x.Site'='ID'))

# now merge the surrogate gauge for the monitoring set to the raw TP df:

df.CSI_TP<-left_join(df.CSI_TP, df.dist[,c(10,9)], by = 'Monitoring.Set')

# then merge the flows based on the surrogate gauge and date:

df.CSI_TP<-left_join(df.CSI_TP, df.CSI.Surro_gauges_Q[,2:4], by = c('site_no', 'Date'))

# some of the monitoring sets remain that didnt have sites with over 50 samples
# filtering these out:

df.CSI_TP<-df.CSI_TP%>%filter(Monitoring.Set %in% df.dist$Monitoring.Set)

# some of the flows are NA. These should be where the flow record didnt overlap the sample record
# so I need to estimate flows on these days for these gauges using another surrogate gauge!
# I think Fall Creek is the one to use:
# I already determined that OLS and DA scaling produce similar results for scaling flow from Fall Creek to the other three surrogate gauges
# so I can use Fall Creek and the DA ratios of each site and Fall Creek to fill in the dates that dont have flow in df.CSI_TP

# determine the drainage areas for the 4 surrogate gauges:

surro_DA<-lapply(unique(df.dist$site_no),\(i) readNWISsite(i)$drain_area_va)%>%
  purrr::set_names(unique(df.dist$site_no))%>%
  bind_rows(., .id = 'site_no')%>%
  pivot_longer(cols = everything(), names_to = 'site_no', values_to = 'Surro_DA_sqmi')

# calculate the drainage area ratios:

surro_DA$Surro_DA_ratio_with_Fall_Creek<-surro_DA$Surro_DA_sqmi/surro_DA$Surro_DA_sqmi[2] # the second one is fall creek

# add these to the TP df:

df.CSI_TP<-left_join(df.CSI_TP,surro_DA, by = "site_no")

# now add the flows at fall creek for every observation:
# first create a df of just the fall creek flows:

fc_Q<-df.CSI.Surro_gauges_Q%>%filter(site_no == "04234000")%>%
  rename(FC_flow = 4)

# now add this as a column to TP df, only use Date to merge:

df.CSI_TP<-left_join(df.CSI_TP, fc_Q[,3:4], by = 'Date')

# then fill in the missing flows using conditional to test if its NA then make the cell value equal to the product of the FC flow column and the DA ratio:

df.CSI_TP<-df.CSI_TP%>%
  mutate(X_00060_00003 = ifelse(is.na(X_00060_00003), Surro_DA_ratio_with_Fall_Creek*FC_flow, X_00060_00003))

# looks like it worked!

# now scale these surrogate gauge flows using the DA of the site. To do this:

# delinate the the CSI sites: I dont need all ofthem right now but might as well do all of them

l.SS_WS.CSI_TP<-lapply(seq_along(df.CSI_TP_sites$ID), \(i) Ryan_delineate(df.CSI_TP_sites$Longitude[i], df.CSI_TP_sites$Latitude[i]))

l.sp.SS_WS.CSI_TP<-l.SS_WS.CSI_TP

for (i in seq_along(l.sp.SS_WS.CSI_TP)){
  
  tryCatch({
    
    l.sp.SS_WS.CSI_TP[[i]]<-toSp(watershed = l.sp.SS_WS.CSI_TP[[i]], what = 'boundary')
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# set the names of the list:

names(l.sp.SS_WS.CSI_TP)<-df.CSI_TP_sites$ID

# need to remove the drainage areas that did not work in toSp (nothing we can do about losing these, idk why some dont come out of the function right):

l.sp.SS_WS.CSI_TP<-l.sp.SS_WS.CSI_TP[sapply(l.sp.SS_WS.CSI_TP, function(x) class(x) == "SpatialPolygonsDataFrame")]

# convert from sp (old gis in r) to sf (new gis in r)

l.sf.SS_WS.CSI_TP<-lapply(l.sp.SS_WS.CSI_TP, st_as_sf)

# check validity of sf objects and then make valid:

lapply(l.sf.SS_WS.CSI_TP, st_is_valid)

l.sf.SS_WS.CSI_TP<-lapply(l.sf.SS_WS.CSI_TP, st_make_valid)

# need to convert the Shape_Area column to numeric for all dfsin the list or bind_rows wont work:

l.sf.SS_WS.CSI_TP<-lapply(l.sf.SS_WS.CSI_TP, \(i) i%>%mutate(Shape_Area = as.numeric(Shape_Area)))

# create a single sf df with all the sample site draiange areas:

df.sf.CSI_TP<-bind_rows(l.sf.SS_WS.CSI_TP, .id = 'Name')%>%
  relocate(Name, .before= 1)

# convert to vect and add a drainage area in km2 column

vect.CSI_TP<-vect(df.sf.CSI_TP)

vect.CSI_TP$area_KM2<-expanse(vect.CSI_TP, unit="km")

# now filter out the downstream sites into a seperate df:

df.sf.CSI_TP_7<-df.sf.CSI_TP%>%filter(Name%in%keep)

mapview(df.sf.CSI_TP_7, zcol = 'Name')+mapview(map.CSI.TP_sites)




# now add the surrogate gauges for the Seneca-Keuka sites:
# I looked at the SK sites in the code SK_CQ.R to prep for comparing SWAT and actual CQ curves.
# Both DOW and CSI were quiered for TPsites in the SK basin, but Iwas using a threshhold of 7samples
# also, with a threshold of 7 samples, there are only 6 subwatersheds with 'good data' meaning that the other monitoring locations are inline with the two lake tributaries (not what we want)
# So going through that code and changing it to a 20sample threshold, we lose two of the watersheds. 
# the 4 that remain are:

SK_keep<-c("Big Stream - Big Stream Mouth @ Glenora Point",
           "Kashong Creek - Bridge at Route 14",
           "Reeder Creek - Reeder Creek Mouth",
           "Catharine Creek - Upstream of Montour Falls WWTP")

# and the surrogate gauges for these I also determined in that other code to be:

# Cathrine Creek will use its gauge
# Big Stream will use Cathrine Creek
# Kashong Creek and Reeder Creek will use Sugar Creek

# however, these surrigate gauges needed to be extended using a long term gauge, Flint Creek (adjacent to SK watershed)

# in the previous code mentioned above, this was done using OLS
# and the resulting C-Q dataframe assembled. That df is read in here:

load("C:/PhD/Research/Data/Seneca-Keuka/df.SK.TP_CQ.Rdata")

# extract lat and long columns from the geometry column
# return a distinct dataframe on the sites:

df.SK.TP_CQ_sites<-tidyr::extract(df.SK.TP_CQ, geometry, into = c('Long', 'Lat'), '\\((.*),(.*)\\)', conv = T)%>%
  distinct(Name, .keep_all = T)

# The watersheds we delinated in the last code, but looking at the Streamstats sf dataframe that was saved,
# it is clear the watershed attributes didn't download (idk, sometimnes SS just returns zeros for all the watershed parameters)
# So I am going to re-delinate these watersheds here to see if I can get those parameters. 
# I will do the delination all at once when I delinate all of the CSI sites.

