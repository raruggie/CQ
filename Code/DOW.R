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

df.NWIS.Q_sites<-read.csv("C:/PhD/CQ/Raw_Data/df.NWIS.Q_sites.csv", colClasses = c(site_no = "character"))

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

# find the sites within X miles using filter on the distance column. Also make the usgs site_no character:

df.dist<-df.dist%>%
  group_by(y.Site)%>%
  filter(dist_meters < 50/meters_to_miles)%>%
  ungroup()%>%
  mutate(x.Site = as.character(x.Site))

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

df.DOW_TP.Surro<-left_join(df.DOW_TP.keep, df.dist, by = c('site_id'='y.Site'))

# download the raw flow data for NWIS sites:
# this is major. There are over 2000 sites in NYS with daily flow data
# there are also over 1000 sites within the 50 mile radius for these DOW sites
# This does notseems worth it. 
# instead I am just going to use the rawdaily flow data i have forthe TP sites (it ended up working out to give some result here, meaning Ifound a surrogate gauge within 50 miles. But maybe if I had the raw daily flow data for all of the 2000 plus NWIS sites, I could find much better gauges)

df.NWIS_Q<-read.csv("C:/PhD/CQ/Raw_Data/df.NWIS.Q.csv", colClasses = c(site_no = "character"))

# filter flowdata for the sites in df.dist and convert date to date (takes a minute with a 50 mile radius). Also change site_no to character:

df.NWIS_Q.DOW<-df.NWIS_Q%>%
  filter(site_no %in% unique(df.dist$x.Site))%>%
  mutate(Date = as.Date(Date), site_no = as.character(site_no))

# merge this dataframe with the DOW TP sample dataframe, which is possible because the potential surrogate gauges were added a few steps before: 

df.DOW_TP.Surro<-left_join(df.DOW_TP.Surro, df.NWIS_Q.DOW, by = c('x.Site'='site_no', 'sample_date'='Date'))

# now remove the NA in the added NWIS columns, which are where there was no flow data for the DOW sample, 
# group by the DOW site and NWIS site (combinaiton of the two) and summarize to get the number paired CQ observations that would existif using that surrogate gauge for that DOW site
# create a column using this n that is the ratio of the totalN for the site to the paired CQ n:

df.DOW_TP.Surro.summarized<-df.DOW_TP.Surro%>%
  drop_na(X_00060_00003)%>%
  group_by(site_id, x.Site)%>%
  summarise(Paired_CQ_n=n())%>%
  arrange(site_id, desc(Paired_CQ_n))%>%
  ungroup()%>%
  left_join(.,df.DOW_TP_sites[,c(1,2)], by = 'site_id')%>%
  mutate(Paired_CQ_n = Paired_CQ_n/n)

# now plot this: to dothis:

# create a mapping df for the NWIS sites in the above df:

map.DOW_TP.Surro.summarized<-df.DOW_TP.Surro.summarized%>%
  select(c(2,3))%>%
  left_join(.,distinct(df.dist[,c(1:3)]), by = 'x.Site')%>%
  st_as_sf(.,coords=c('x.Longitude','x.Latitude'), crs = 4326)

# make a map (already have the DOW TPsample sites): to do this:

# add 100 mile radius to these DOW points to see which sites are for it:

dat_circles <- st_buffer(map.DOW_TP.Surro.summarized, dist = 50/meters_to_miles)

# now map:

mapview(map.DOW.TP.sites_keep, col.regions=list("darkred"))+mapview(dat_circles,col.regions=list("red"))+ mapview(map.DOW_TP.Surro.summarized, zcol= 'Paired_CQ_n')

# it looks like at 50 miles, each DOW site has a gauge with an actual paired to potential total paired CQ ratio of 1
# however, it is possible that a yellow (ratio = 1) in the plot above is not for the site even if it is in its radius since radi can overlap
# but we cancheck to make sure that each DOW site is represented after filtering the df to that with a ratio of 1:

df.DOW_TP.Surro.summarized.final<-df.DOW_TP.Surro.summarized%>%filter(Paired_CQ_n==1)

unique(df.DOW_TP.Surro.summarized.final$site_id)

# all 9 sites are represented, mean each one has a surogate gauge with a ratio of 1!

# and now relook at the map but in the pltting df keep the DOW site_id to align them:

# create map variable:

map.DOW_TP.Surro.summarized.final<-df.DOW_TP.Surro.summarized.final%>%
  select(c(1,2))%>%
  left_join(.,distinct(df.dist[,c(1:3)]), by = 'x.Site')%>%
  st_as_sf(.,coords=c('x.Longitude','x.Latitude'), crs = 4326)

# put map together:

mapview(map.DOW.TP.sites_keep, zcol = 'site_id')+
  mapview(dat_circles, col.regions=list("red"), alpha.regions = 0.2)+
  mapview(map.DOW_TP.Surro.summarized.final, zcol= 'site_id')

# from this map the best USGS gauges for the sites in the order are:

v.DOW.TP.sites_surro<-sort(c("02-ALGY-20.3:04213500", '06-SUSQ-6.9:01531000','14-DELA-1.3:01421000',"06-NVUS-0.9:01421618", "12-MOHK-136.0:04250200", "12-ORSK-0.9:01349150", "11-UHUD-98.3:01327750", "11-UHUD-64.0:01327750", "09-CGAY-2.7:04269000"))

# transform this vector in a df:

df.DOW.TP.sites_surro<-data.frame(do.call(rbind, strsplit(v.DOW.TP.sites_surro, ":", fixed=TRUE)))%>%
  rename(site_id = 1, site_no_keep = 2)

# filter down df.DOW_TP.keep to just the CQ paired observations. To do this:

# merge the chosen surrogate gauge name to each DOW observation, then filter where the gauge observations are:

df.DOW_TP.keep.final<-df.DOW_TP.Surro%>%
  left_join(., df.DOW.TP.sites_surro, by = 'site_id')%>%
  filter(x.Site == site_no_keep)

# now add in columns to draiange area scale surrogate flows: to do this:

# create a df of the unique DOW sites

temp<-df.DOW_TP_sites%>%
  filter(site_id %in% v.DOW.TP.sites_keep)

# delinate the DOW sites and convert to sf unsing functions (sourced):

l.SS_WS.DOW<-lapply(seq_along(temp$site_id), \(i) Delineate(temp$longitude[i], temp$latitude[i]))

df.sf.DOW.DA<-fun.l.SS_WS.to.sfdf(l.SS_WS.DOW,temp$site_id)

# look at a map of the watersheds:

mapview(df.sf.DOW.DA, zcol = 'Name')

# convert to vect and add a drainage area in km2 column

vect.DOW<-vect(df.sf.DOW.DA)

vect.DOW$area_KM2<-expanse(vect.DOW, unit="km")

# add DA in sqmi to raw TP data: to do this:

# add it to temp so that you have unique sample location to its DA:

temp$DA_sqmi<-vect.DOW$area_KM2*0.386102

# now merge with raw TP data:

df.DOW_TP.keep.final<-left_join(df.DOW_TP.keep.final, temp[,c(1,6)], by = 'site_id')

# now get the surrogate drainage areas:

temp<-readNWISsite(unique(df.DOW_TP.keep.final$x.Site))[,c(2,30)]

# merge the DA with the raw TP data, calcuate the DA ratio, and scale the flows:

df.DOW_TP.keep.final<-left_join(df.DOW_TP.keep.final, temp, by = c("x.Site"="site_no"))%>%
  mutate(DA_ratio=DA_sqmi/drain_area_va)%>%
  mutate(Q_scaled = X_00060_00003*DA_ratio)

# now I have CQ curves and SS metadata

# plots of CQ curves:
  
ggplot(df.DOW_TP.keep.final, aes(x = log(Q_scaled), y = log(result_value)))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap(dplyr::vars(site_id), scales = 'free')

# because I dont trust the delinations, I am going to wait to use them to get land use, elevation, climate predictors 

# save.image(file = 'C:/PhD/CQ/Processed_Data/DOW.Rdata')

