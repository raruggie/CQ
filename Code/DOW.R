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

# combine the two dataframes:

map.NWIS_and_DOW.TP_sites<-bind_rows(map.DOW.TP_sites,map.NWIS.TP_sites)

# then map:

mapview(map.NWIS_and_DOW.TP_sites, zcol = 'agency_cd')

# looking at this map, maybe there are a hand full of potential useful sites in the DOW database
# to include in the CQanalysis, but for the most part  I think the NWIS sites cover the same watersheds

# I did pair up flow data when the sample threshold was 99. But at 20 I dont want to do it
# since the workflow is by hand.

# below is the workflow I used when the threshold was at 99 samples:

# filter the NYS NWIS sites to those with over 1000 flow days:

df.NWIS<-df.NWIS.Q%>%ungroup()%>%mutate(site_id = site_no, n = nflowdays)%>%
  select(site_id, n,latitude, longitude, agency_cd)%>%
  filter(n>1000)

# merge the TP sites with its lat long metadata
# also filter to sites with over 99 samples and bind rows with the NWIS sites 
# and then plot:

x<-left_join(df.temp,sms[,c(1,5,6)], by = 'site_id')%>%
  mutate(agency_cd = 'DOW Sample')%>%
  filter(n>20)%>%
  # bind_rows(df.NWIS)%>%
  drop_na(latitude,longitude)%>%
  st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)%>%
  mapview(., zcol = 'agency_cd')

# the following DOW-USUS site pairs:

# 07-OSWE-5.2 04249000
# 
# 08-BLCK-1.4 04260500 
# 
# 11-UHUD-64.0 01326500
# 
# 12-MOHK-1.5 	01357500 # this site looks potential for backwater from adjacent Husdon river to impact it.
# 
# 11-UHUD-2.7 01335754
# 
# 13-LHUD-120.2 01358000 # this gauge is quite far upstream from sample
# 
# 13-LHUD-66.3  01358000 # there isnt a gauge anywhere near here. Husdon river at green island is the last gauge.
# 
# 14-DELA-1.3 	01428500 # the gaue is also on the delware river butthe monguap flows into it before the sample locations IO think it might just be betterto use CBNTN sites if we aregoing to use delware data.
# 
# 06-NANG-0.7 01512500 # this sample site is right in binghamton urban area
# 
# 06-SUSQ-6.9 01513831
# 
# 01-BUFF-1.7 04214500 # caz creek flows in before sample site
# 
# 04-GENS-2.6 04232000 # proably would just use USGS sample data for gene river
# 
# 07-SEOS-22.4 04237496

# do any of these gauge sites show up in the NWIS analysis from above?

DOW_USGS_surrogates<-c('04249000','04260500','01326500','01357500','01335754','01358000','01358000','01428500','01512500','01513831','04214500','04232000','04237496')

DOW_USGS_surrogates<-df.NWIS.TP_CQ%>%filter(site_no %in% DOW_USGS_surrogates)%>%
  distinct(site_no)

# four of these sites pop up.
# lets look at them on the map:

df.NWIS<-df.NWIS%>%filter(site_id %in% DOW_USGS_surrogates$site_no)

left_join(df.temp,sms[,c(1,5,6)], by = 'site_id')%>%
  mutate(agency_cd = 'DOW Sample')%>%
  filter(n>50)%>%
  bind_rows(df.NWIS)%>%
  drop_na(latitude,longitude)%>%
  st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)%>%
  mapview(., zcol = 'agency_cd')

# these are pretty much all colocate, and the NWIS has more data
# these will be filtered these out

char_vector<-c('07-OSWE-5.2', '04249000','08-BLCK-1.4','04260500',
               '11-UHUD-64.0','01326500','12-MOHK-1.5','01357500',
               '11-UHUD-2.7','01335754','13-LHUD-120.2','01358000',
               '13-LHUD-66.3',  '01358000',
               '14-DELA-1.3', 	'01428500' ,
               '06-NANG-0.7', '01512500' ,
               '06-SUSQ-6.9', '01513831',
               '01-BUFF-1.7', '04214500' ,
               '04-GENS-2.6', '04232000' ,
               '07-SEOS-22.4', '04237496')

odd_elements <- char_vector[(1:length(char_vector)) %% 2 == 1]
even_elements <- char_vector[(1:length(char_vector)) %% 2 == 0]

# Create a dataframe and filter the (nearly) colocated sites:

df.DOW_TP_sites_metadata<-data.frame(DOW = odd_elements, Surrogate_Gauge = even_elements)%>%
  filter(!Surrogate_Gauge %in% DOW_USGS_surrogates$site_no)

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
