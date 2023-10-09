# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

library(climateR)
library(geosphere)
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

# geosphere functions return meters. If you want miles:

meters_to_miles = 1/1609.334

####################### Functions #######################

source("C:/PhD/CQ/Code/Ryan_functions.R")

####################### Goal of code #######################

# Process the DOW database for the CQ analysis

####################### Workflow #######################

# the DOW portal is interesting because we have lots of different
# consituents measured. But the issue is that we wont necessiary have
# overlapping observations, and we wont necessairly have these consituents measured 
# at other sites, as to build a RF model to predict P using many other consitents in the sample.

# import in the DOW raw sample and site metadata csv files

sc<-read.csv('C:/PhD/CQ/Raw_Data/DOW_Streams_Chemistry.csv')
sms<-read.csv('C:/PhD/CQ/Raw_Data/DOW_Stream_Monitoring_Sites.csv') 

# filter chemistry for TP:

df.DOW_TP<-sc%>%filter(parameter_name=="PHOSPHORUS, TOTAL")

# I want to look at how many samples are there per unique site:

df.DOW_TP_sites<-df.DOW_TP%>%group_by(site_id)%>%summarize(n=n())%>%
  arrange(desc(n))

# I want to see if these sites match up with a USGS gauge, but this is a pain in the ass!
# some sites might fall on the same stream that has gauge data, which would be easy to then use DA scaling to get the flows at the sample site
# but the workflow to do this is to do it by HAND (unless I want to automate it, which I could proablyfigure out but it wouldnt be perfect and would take time to learn how to do it)
# also, if the sample location is not on a stream, I would need to come up with a suitable nearby gauge, which again would be by HAND

# lets just start with making a df of DOW sites with number of TP samples > 20:

df.DOW_TP_sites<-left_join(df.DOW_TP_sites,sms[,c(1,5,6)], by = 'site_id')%>%
  mutate(agency_cd = 'DOW Sample')%>%
  filter(n>20)

# and look at the map:

map.DOW.TP_sites<-df.DOW_TP_sites%>%
  drop_na(latitude,longitude)%>%
  st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)

mapview(map.DOW.TP_sites, zcol = 'agency_cd')

# lets add to this map the NWIS TP sites: to do this

# import the NWIS site map:

load('C:/PhD/CQ/Processed_Data/map.NWIS.TP_sites.Rdata')

# combine the two dataframes:

map.NWIS_and_DOW.TP_sites<-bind_rows(map.DOW.TP_sites,map.NWIS.TP_sites)

# then map:

mapview(map.NWIS_and_DOW.TP_sites, zcol = 'agency_cd')

# looking at this map, maybe there are a hand full of potential useful sites in the DOW database
# to include in the CQanalysis, but for the most part  I think the NWIS sites cover the same watersheds

# from this map, the following DOW sites stand out as different from the NWIS sites
# and useful, meaning not in an urban setting:

v.DOW.TP.sites_keep<-sort(c('02-ALGY-20.3','06-SUSQ-6.9', '12-ORSK-0.9', '12-MOHK-136.0', '11-UHUD-98.3', '11-UHUD-64.0', '06-NVUS-0.9', '09-CGAY-2.7', '14-DELA-1.3'))

# now looking at a map of just these sites:

map.DOW.TP.sites_keep<-df.DOW_TP_sites%>%
  filter(site_id %in% v.DOW.TP.sites_keep)%>%
  drop_na(latitude,longitude)%>%
  st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)

mapview(map.DOW.TP.sites_keep, zcol = 'agency_cd')

# I want USGS gauges within a radius of X miles of these sites from the master list of NYS gauges. To do this:

# read in the NYS USGS gauges with daily flow:

df.NWIS.Q_sites<-read.csv("C:/PhD/CQ/Raw_Data/df.NWIS.Q_sites.csv")

# determine all possible combinaitons of the USGS gauges and DOW locations: to do this 

# establish the DOW keep dataframe:

df.DOW.TP.sites_keep<-df.DOW_TP_sites%>%
  filter(site_id %in% v.DOW.TP.sites_keep)

# expand() is used here to create all posible combinaitons of the USGS gauges and the 7 CSI downstream sites:
# then left join the lat longs of the CSI sites (already have the gauge lat longs)

df.dist<- df.NWIS.Q_sites%>%
  group_by(site_no, dec_lat_va, dec_long_va)%>%
  expand(site_id = df.DOW.TP.sites_keep$site_id)%>%
  left_join(.,df.DOW.TP.sites_keep, by = 'site_id')%>%
  dplyr::select(c(1:4,6,7))%>%
  arrange(site_id)

# rename the columns:

names(df.dist)<-c("x.Site","x.Latitude","x.Longitude","y.Site","y.Latitude","y.Longitude")

# add a distance column (distance between every possible combinaiton of USGS gauge and CSI site) by calling distHaversine (vectorized) on each pair:

df.dist$dist_meters <- geosphere::distHaversine(df.dist[3:2], df.dist[6:5]) # units of meters

# find the sites within X miles using filter on the distance column:

df.dist<-df.dist%>%
  group_by(y.Site)%>%
  filter(dist_meters < 50/meters_to_miles)%>%
  ungroup()

# create a df for mapping: to do this:

# add ID columns:

df.dist_for_map<-df.dist%>%
  mutate(x.Type = 'USGS', .after = 3)%>%
  mutate(y.Type = 'DOW', .after = 7)

# remove the .x and .y from the col names (so rbind works)
names(df.dist_for_map) <- substring(names(df.dist_for_map), 3)

# use rbind to stack the USGS and OW sites and lat longs and type columns into one dataframe, then convert to sf df:

map.DOW.Surro_gauges<-rbind(df.dist_for_map[1:4], df.dist_for_map[5:8])%>%
  st_as_sf(.,coords=c('Longitude','Latitude'), crs = 4326)

# map:

mapview(map.DOW.Surro_gauges, zcol = 'Type')

# going through by hand to find the best gauge for each site:
# this is going to be hard because I dont know which gauges have the most flow data.

# I can pair the flow data to the TP sample data and see which sites have the most numberof paired CQ observations:

# filter the raw DOW TP data to the sites we are keeping and convert dates for later merging:

df.DOW_TP.keep<-df.DOW_TP%>%filter(site_id %in% v.DOW.TP.sites_keep)%>%
  mutate(sample_date_time = as.POSIXct(sample_date, format = '%m/%d/%Y, %I:%M %p'), .before = 2)%>%
  mutate(sample_date = as.Date(sample_date_time))

# merge to this data frame the df.dist results (with potential surrogate NWIS gauges paired to DOW sample sites):

df.DOW_TP.keep<-left_join(df.DOW_TP.keep, df.dist, by = c('site_id'='y.Site'))

# read in the raw flow data for NWIS sites:

# df.NWIS_Q<-read.csv("C:/PhD/CQ/Raw_Data/df.NWIS.Q.csv")

# filter for the sites in df.dist and convert date to date:

df.NWIS_Q.DOW<-df.NWIS_Q%>%
  filter(site_no %in% unique(df.dist$x.Site))%>%
  mutate(Date = as.Date(Date))

# merge this dataframe with the DOW TP sample dataframe, which is possible because the potential surrogate gauges were added a few steps before: 

df.DOW_TP.keep<-left_join(df.DOW_TP.keep, df.NWIS_Q.DOW, by = c('x.Site'='site_no', 'sample_date'='Date'))

# now remove the NA in the added NWIS columns, which are where there was no flow data for the DOW sample, 
# group by the DOW site and NWIS site (combinaiton of the two) and summarize to get the number paired CQ observations that would existif using that surrogate gauge for that DOW site
# create a column using this n that is the ratio of the totalN for the site to the paired CQ n:

df.DOW_TP.keep.summarized<-df.DOW_TP.keep%>%
  drop_na(X_00060_00003)%>%
  group_by(site_id, x.Site)%>%
  summarise(Paired_CQ_n=n())%>%
  arrange(site_id, desc(Paired_CQ_n))%>%
  ungroup()%>%
  left_join(.,df.DOW_TP_sites[,c(1,2)], by = 'site_id')%>%
  mutate(Paired_CQ_n = Paired_CQ_n/n)

# now plot this: to dothis:

# create a mapping df for the NWIS sites in the above df:

map.DOW_TP.keep.summarized<-df.DOW_TP.keep.summarized%>%
  select(c(2,3))%>%
  left_join(.,distinct(df.dist[,c(1:3)]), by = 'x.Site')%>%
  st_as_sf(.,coords=c('x.Longitude','x.Latitude'), crs = 4326)

# make a map (already have the DOW TPsample sites): to do this:

# add 100 mile radius to these DOW points to see which sites are for it:

dat_circles <- st_buffer(map.DOW.TP.sites_keep, dist = 50/meters_to_miles)

# now map:

mapview(map.DOW.TP.sites_keep, col.regions=list("darkred"))+mapview(dat_circles,col.regions=list("red"))+ mapview(map.DOW_TP.keep.summarized, zcol= 'Paired_CQ_n')

# it looks like at 50 miles, each DOW site has a gauge with an actualpaired to potential total paired CQ ratio of 1
# however, it is possible that a yellow (ratio = 1) in the plot above is not for the site even if it is in its radius since radi can overlap
# but we cancheck to make sure that each DOW site is represented after filtering the df to that with a ratio of 1:

df.DOW_TP.keep.summarized.final<-df.DOW_TP.keep.summarized%>%filter(Paired_CQ_n==1)

unique(df.DOW_TP.keep.summarized.final$site_id)

# all 9 sites are represented, mean each one has a surogate gauge with a ratio of 1!

# and now relook at the map but in the pltting df keep the DOW site_id to align them:

map.DOW_TP.keep.summarized.final<-df.DOW_TP.keep.summarized.final%>%
  select(c(1,2))%>%
  left_join(.,distinct(df.dist[,c(1:3)]), by = 'x.Site')%>%
  st_as_sf(.,coords=c('x.Longitude','x.Latitude'), crs = 4326)

mapview(map.DOW.TP.sites_keep, zcol = 'site_id')+
  mapview(dat_circles, col.regions=list("red"), alpha.regions = 0.2)+
  mapview(map.DOW_TP.keep.summarized.final, zcol= 'site_id')

# from this map the best USGS gauges for the sites in the order are:

v.DOW.TP.sites_surro<-sort(c("02-ALGY-20.3:4213500", '06-SUSQ-6.9:1531000','14-DELA-1.3:1421000',"06-NVUS-0.9:1421618", "12-MOHK-136.0:4250200", "12-ORSK-0.9:1349150", "11-UHUD-98.3:1327750", "11-UHUD-64.0:1327750", "09-CGAY-2.7:4269000"))

# transform this vector in a df:

df.DOW.TP.sites_surro<-data.frame(do.call(rbind, strsplit(v.DOW.TP.sites_surro, ":", fixed=TRUE)))%>%
  rename(site_id = 1, site_no = 2)

# the good thing now is that I have the paired CQ data!


save.image(file = 'C:/PhD/CQ/Processed_Data/DOW.Rdata')







# now download the raw flow data for these gauges:

df.DOW_surro_Qs<-readNWISdv(siteNumbers = df.DOW_TP_sites_metadata$Surrogate_Gauge, parameterCd = '00060')

# then download the metadata as to get their drainage areas

df.DOW_surro_Qs_metadata<-readNWISsite(siteNumbers = df.DOW_TP_sites_metadata$Surrogate_Gauge)

# then left join with the DOW site name and convert from mi to km2:

df.DOW_TP_sites_metadata<-left_join(df.DOW_TP_sites_metadata, df.DOW_surro_Qs_metadata[,c(2,30)], by = c('Surrogate_Gauge'='site_no'))%>%
  rename(Surrogate_DA = 3)%>%
  mutate(Surrogate_DA=Surrogate_DA*2.58999)

# then delinate the drainage areas of the sampling sites to get their drianage areas:

# first add the lat longs to the metadata df:

df.DOW_TP_sites_metadata<-left_join(df.DOW_TP_sites_metadata,sms[,c(1,5,6)],by = c('DOW'='site_id') )

# next delinate using lapply to get into list

l.sf.SS_WS.DOW_TP_sites[[2]]<-x

l.SS_WS.DOW_TP_sites<-lapply(seq_along(df.DOW_TP_sites_metadata$latitude), \(i) Ryan_delineate(df.DOW_TP_sites_metadata$longitude[i],df.DOW_TP_sites_metadata$latitude[i]))

# then convert out of watershed to spatial polygons df: need to use a loop because toSp function doesnt work always with lapply

l.sp.SS_WS.DOW_TP_sites<-l.SS_WS.DOW_TP_sites

for (i in seq_along(l.sp.SS_WS.DOW_TP_sites)){
  
  tryCatch({
    
    l.sp.SS_WS.DOW_TP_sites[[i]]<-toSp(watershed = l.sp.SS_WS.DOW_TP_sites[[i]], what = 'boundary')
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# first set the names of the list:

names(l.sp.SS_WS.DOW_TP_sites)<-df.DOW_TP_sites_metadata$DOW

# need to remove the drainage areas that did not work in toSp (nothing we can do about losing these, idk why some dont come out of the function right):

l.sp.SS_WS.DOW_TP_sites<-l.sp.SS_WS.DOW_TP_sites[!sapply(l.sp.SS_WS.DOW_TP_sites, function(x) class(x) == "watershed")]

# convert from sp (old gis in r) to sf (new gis in r)

l.sf.SS_WS.DOW_TP_sites<-lapply(l.sp.SS_WS.DOW_TP_sites, st_as_sf)

# check validity of sf objects and then make valid:

lapply(l.sf.SS_WS.DOW_TP_sites, st_is_valid)

l.sf.SS_WS.DOW_TP_sites<-lapply(l.sf.SS_WS.DOW_TP_sites, st_make_valid)

# need to convert the Shape_Area column to numeric for all dfsin the list or bind_rows wont work:

l.sf.SS_WS.DOW_TP_sites<-lapply(l.sf.SS_WS.DOW_TP_sites, \(i) i%>%mutate(Shape_Area = as.numeric(Shape_Area)))

# create a single sf df with all the sample site draiange areas:

df.sf.DOW_TP_sites<-bind_rows(l.sf.SS_WS.DOW_TP_sites, .id = 'Name')%>%
  relocate(Name, .before= 1)

# look at a map of the watersheds:

mapview(df.sf.DOW_TP_sites, zcol = 'Name')

# the two hudon river ones didnt work
# I tried delinating them a second time but still gives zeros in the watershed attributes columns

# convert to vect and add a drainage area in km2 column

vect.DOW_TP_sites<-vect(df.sf.DOW_TP_sites)

vect.DOW_TP_sites$area_KM2<-expanse(vect.DOW_TP_sites, unit="km")

# create a temp dataframe of the DOW site names and the draiange areas out of the SpatVect object:

temp<-as.data.frame(vect.DOW_TP_sites[,c('Name', 'area_KM2')])

# now merge these drainage area calcs with metadata df:
# and calcuate the DA ratios:

df.DOW_TP_sites_metadata<-left_join(df.DOW_TP_sites_metadata, temp, by = c('DOW'='Name'))%>%
  rename(DOW_DA = 6)%>%
  mutate(DA_ratio = DOW_DA/Surrogate_DA)

df.DOW_TP_sites_metadata<-df.DOW_TP_sites_metadata%>%
  mutate(DA_ratio = DOW_DA/Surrogate_DA)

# I want to look at the drainage areas of the gauges these to make sure they are right:

temp<-lapply(seq_along(df.DOW_surro_Qs_metadata$dec_lat_va), \(i) Ryan_delineate(df.DOW_surro_Qs_metadata$dec_long_va[i], df.DOW_surro_Qs_metadata$dec_lat_va[i]))

for (i in seq_along(temp)){
  
  tryCatch({
    
    temp[[i]]<-toSp(watershed = temp[[i]], what = 'boundary')
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# first set the names of the list:

names(temp)<-df.DOW_TP_sites_metadata$Surrogate_Gauge

# need to remove the drainage areas that did not work in toSp (nothing we can do about losing these, idk why some dont come out of the function right):

temp<-temp[!sapply(temp, function(x) class(x) == "watershed")]

# convert from sp (old gis in r) to sf (new gis in r)

temp<-lapply(temp, st_as_sf)

# check validity of sf objects and then make valid:

lapply(temp, st_is_valid)

temp<-lapply(temp, st_make_valid)

# need to convert the Shape_Area column to numeric for all dfsin the list or bind_rows wont work:

temp<-lapply(temp, \(i) i%>%mutate(Shape_Area = as.numeric(Shape_Area)))

# create a single sf df with all the sample site draiange areas:
temp<-bind_rows(temp, .id = 'Name')%>%
  relocate(Name, .before= 1)

# look at a map of the watersheds:

mapview(df.sf.DOW_TP_sites, zcol = 'Name')+mapview(temp, zcol = 'Name')

# they look fine (fyi some didnt show up, but the others check out so i think its ok)

# next pair up the flows to TP observations and scale the flows:

# Steps to do this:
# filter the DOW_TP obsevations to the sites detemrined above (have over 50 samples)
# convert dattime to POSIXct and create a column for just the date
# left join to add the surrogate gauge data and DA ratio
# left join to add the flow data by site and date
# rename the surrogate gauge flow column name and calcute the scaled flows by multiplying by the DA_ratio
# remove sites/observations where flow and TP didnt fall on same date:

df.DOW_TP_CQ<-df.DOW_TP%>%
  filter(site_id %in% df.DOW_TP_sites_metadata$DOW)%>%
  mutate(sample_date = as.POSIXct(sample_date, format = '%m/%d/%Y, %I:%M %p'))%>%
  mutate(Date = as.Date(sample_date), .after = 2)%>%
  left_join(., df.DOW_TP_sites_metadata, by = c('site_id'='DOW'))%>%
  left_join(., df.DOW_surro_Qs[,c(2:4)], by = c('Surrogate_Gauge'='site_no', 'Date'='Date'))%>%
  rename(Surrogate_Gauge_Flow = X_00060_00003)%>%
  mutate(Q_DA_scaled = Surrogate_Gauge_Flow*DA_ratio)%>%
  drop_na(Q_DA_scaled)%>%
  filter(Q_DA_scaled>0)

# now make plots of CQ curves:

ggplot(df.DOW_TP_CQ, aes(x = log(Q_DA_scaled), y = log(result_value)))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap(dplyr::vars(site_id), scales = 'free')

# Done! The df df.DOW_TP_sites_metadata serves as the CQ and metadata dfs
# since it has the lat longs for each site
