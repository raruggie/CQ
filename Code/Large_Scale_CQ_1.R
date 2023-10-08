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

source("C:/PhD/Research/Code/Manuscript/temp/Ryan_functions.R")

####################### Goal of code #######################

# this is the first code of two files to set up a dataframe 
# for multi-site, multi-basin CQ analysis 

# in this code, the lat/longs, C, and Q data are assembled

# in the next file, the watershed attributes will be calcualted from the
# lat/longs and datetimes of CQ paired observations

# Steps for this code:

####################### Workflow #######################

####------------Step 1 Determine which sites to include------------####

# there are a few main repositories to get data from that I have explored so far:
# NWIS, DOW, CSI, Charley Driscoll

# there is also CEAP but I havnt wored with that data yet
# there is also databases in the papers I have read that we could scour for more data

# for now lets stick with the ones I have worked with, NWIS, DOW, CSI, and Charley (Skan)

# NWIS: is a statewide (potentially region wide) repositiory 
# that is easy to pair flow data with.
  # I already queried this database for sites with over 99 samples and daily flow data

# DOW: unless colocated, need to pair up surrogate flow site.

# CSI: unless colocated, need to pair up surrogate flow site. 

# Skan: Already decided that Owasco inlet is the surrogate flow gauge
  
####------------Step 2: NWIS Import------------####

# this code adapted from NYS_site_ranking.R

#### Part 1: Download the metadata for sites with daily flow data

# Use dataRetrieval::whatNWISdata to query USGS gauge sites with daily flow data in NY?

# USGS_gauge_sites_in_NY

df.NWIS.Q_sites<- whatNWISdata(stateCd = 'NY', parameterCd = "00060")

# clean up this dataframe

# there are duplicates of site numbers with different numbers of daily flow observations
# so group by the site no numbers and just keep the observation with the highest number of flow observations
# I realize now that this is proably not needed since I am grouping by the site number so the resulting number used in the
# raw flow data download will have what it has, but I am just goin to keep this (aka i could just use distinct?)
# also rename some columns and convert them to numeric:

df.NWIS.Q_sites<-df.NWIS.Q_sites%>%
  group_by(site_no)%>%
  slice(which.max(count_nu))%>%
  rename(latitude = dec_lat_va, longitude = dec_long_va, nflowdays = count_nu, begin_date_flow = begin_date, end_date_flow = end_date)%>%
  mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# NOTE: whatNWISsites(stateCd = "NY", parameterCd = "00060") returns a dataframe with rnow = 1465
# while whatNWISdata(stateCd = "NY", parameterCd = "00060") returns a dataframe with nrow = 2249 ?!?!
# this is because there are duplicates. Using the group_by, slice functions gives a dataframe with
# nrow = 1452 for the whatNWISdata thing.

#### Part 2: From the sites with daily flow data, determine which sites also have over 
# 99 discrete samples from one or more of the constituents of interests

# use function I created (sourced) to get dataframe of sites for just TP:

df.NWIS.TP<-fun.df.Pair_consit_flow('00665', df.NWIS.Q_sites)

# download the raw daily flow data for these sites (this takes a looong time):

df.NWIS.Q<-readNWISdv(siteNumbers = df.NWIS.TP$site_no, parameterCd = '00060', startDate = "", endDate = "", statCd = "00003")

# download the raw discrete TP sample data:

df.NWIS.TP<-readNWISqw(siteNumbers = df.NWIS.TP$site_no, parameterCd = '00665')

# now join the TP with the flow df:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP, df.Q, by=c("site_no"="site_no", "sample_dt"="Date"))

# remove observations where there are not CQ pairs:

df.NWIS.TP_CQ<-df.NWIS.TP_CQ%>%drop_na(X_00060_00003)

# arrange by number of TP observations. To do this:

# first arrange the dataframe with the number of samples: 

df.NWIS.TP<-df.NWIS.TP%>%arrange(desc(count_nu))

# then add a column with the rank:

df.NWIS.TP$count_nu_rank<-1:nrow(df.NWIS.TP)

# finally merge this df with df.NWIS.TP_CQ and arrange by the new column:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP_CQ,df.NWIS.TP[,c(1,13)], by='site_no')%>%
  arrange(count_nu_rank)

# next step is to use Q yield to get better looking plot? idk if it will help but want to try
# also doesnt hurt to have the watershed areas as well. To do this:

# download the drainage areas from site metadata using readNWISsite

df.NWIS.TP_site_metadata<-readNWISsite(siteNumbers = unique(df.NWIS.TP_CQ$site_no))

# then select just the site number and DA column:

df.DA<-df.NWIS.TP_site_metadata%>%
  select(site_no, drain_area_va)

# finally merge with df.NWIS.TP_CQ and create a new Q column with area normalized flows (not worrying about units right now): 
# Note: will filter for NA in C and Q for breakpoint analysis in the next step as to keep the full list of sites with CQ pairs in this dataframe.
# Note: Some sites returned NA on draiange areas in readNWISsite, but I'll delinate anyways so I want the full list:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP_CQ, df.DA, by = 'site_no')%>%
  mutate(Q_yield = X_00060_00003/drain_area_va)

# the next step is to fit breakpoint models to the CQ curves. To do this:

# create empty list to hold the results of the segmente function (which does the breakpoint analysis)

l_Seg<-list()

# create a matrix to hold the results of the davies test, which determines if a two slope model is warrented over a single slope model:

davies.test.matrix<-NULL

# create a new dataframe of only paired CQ observations (such that the breakpoit analysis function runs smoothly) (I didnt want to do this in loop for some reason, I cant remeber why but it wouldnt work):

df.NWIS.TP_CQ_for_BP<-df.NWIS.TP_CQ%>%
  drop_na(result_va, Q_yield)

# create a vector of ordered unique site names:

temp<-sort(unique(df.NWIS.TP_CQ_for_BP$site_no))

# test i for for loop building:

# i<-unique(df.NWIS.TP_CQ_for_BP$site_no)[4]

# loop through the sites:

for (i in seq_along(temp)){
  
  tryCatch({
    
    # print the site name for loop debugging:
    
    print(i)
    print(temp[i])
    
    # create a dataframe that will work with segmented. To do this: 
      # filter for the site the loop is in
      # add log transformed C and Q columns, as well as duplicated columns for renamed C and Q
      # filter for real log C and Q values so breakpoint analysis works smoothly:
    
    df<-df.NWIS.TP_CQ_for_BP%>%
      filter(site_no == temp[i])%>%
      mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
      filter(is.finite(log_C))%>%
      filter(is.finite(log_Q))
    
    # build a single slope lm for log C and Q. Tis model is also used inthe breakpoint analysis inthenext step:
    
    m<-lm(log_C~log_Q, df)
    
    # perform breakpoint regression:
    
    m_seg<-segmented(obj = m, npsi = 1)
    
    # perform davies test for constant linear predictor:
    # the results are saved as a string with the the site name and true/false:
    
    x<-paste(temp[i], '-', davies.test(m)$p.val<0.05)
    
    # add the results of davies test to the matrix made prior to this for loop:
    
    davies.test.matrix<-c(davies.test.matrix,x)
    
    # save the breakpoints
    
    bp<-m_seg$psi[1]
    
    # save the slopes: To do this:
    # a conditional statement is needed since sometimes the segmented function wont fit a two slope model at all and will return a object that doesnt work with the slope function used here:
    
    if(length(class(m_seg))==2){
      s<-as.data.frame(slope(m_seg))
    } else{
      s<-NA
    }
    
    # get the intercepts (again conditional statement is needed):
    
    if(length(class(m_seg))==2){
      inter<-as.data.frame(intercept(m_seg))
    } else{
      inter<-NA
    }
    
    
    # get the model fitted data and put in a dataframe:
    
    fit <- data.frame(Q = df$log_Q, Seg_C = fitted(m_seg))
    
    # reformat this dataframe to export out of the loop
    
    if(length(class(m_seg))==2){
      result_df<-fit%>%mutate(site = temp[i], Q_real = df$Q, C = df$C, I1 = inter$Est.[1], I2 = inter$Est.[2], Slope1 = s[1,1], Slope2 = s[2,1], BP = bp)
    } else{
      result_df<-fit%>%mutate(site = temp[i], Q_real = df$Q, C = df$C, I1 = NA, I2 = NA, Slope1 = NA, Slope2 = NA, BP = NA)
    }
    
    l_Seg[[i]]<-result_df
    
  }, 
  
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# Looking at the davies test results:

davies.test.matrix

# transform this matrix into a two column dataframe for use later with df_Seg and plotting. To do this:
# separate the matrix into two columns
# use mutate to remove white space around these character columns:

df.davies<-as.data.frame(davies.test.matrix)%>%
  separate_wider_delim(1, "-", names = c("site", "BP_yes"))%>%
  mutate(across(c(1,2), trimws))

# combine the list of dfs of the breakpoint analysis results (with fitted values, intercepts and slopes) into a single df. To do this:
# merge with the df above to add the BP_yes column
# use replace and the BP yes column with a conditional statement to set the breakpoint Q and C column rows to NA, as to not plot the segmeneted line if davies test was false
# add the number of samples for the site column (for later to plot the sites CQ curves in order):

df_Seg<-bind_rows(l_Seg)%>%
  left_join(., df.davies, by = 'site')%>%
  mutate(across(c(1,2), ~replace(., BP_yes == 'FALSE', NA)))%>%
  left_join(., df.NWIS.TP[,c(1,13)], by = c('site'='site_no'))%>%
  arrange(count_nu_rank)

# ready to plot:

ggplot(df_Seg, aes(x = log(Q_real), y = log(C)))+
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'tomato')+
  facet_wrap(dplyr::vars(count_nu_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Done with NWIS! the two dataframes moving on are 
# df.NWIS.TP_CQ and df.NWIS.TP_site_metadata

# note that this is just for TP!! In the future I may come back to 
# this code and add other consituents

# one last thing: lets look at a map of these:

map.NWIS.TP_sites<-df.NWIS.TP_site_metadata%>%
  rename(longitude=8,latitude=7)%>%
  drop_na(latitude,longitude)%>%
  st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)%>%
  mutate(site_id= paste(site_no, station_nm), n = NA, .before = 1)%>%
  select(c(1:3))

mapview(map.NWIS.TP_sites)

####------------Step 3: DOW Import------------####

# the DOW portal is interesting because we have lots of different
# consituents measured. But the issue is that we wont necessiary have
# overlapping observations, and we wont necessairly have these consituents measured 
# at other sites, as to build a RF model to predict P using many other consitents in the sample.

# import raw DOW database site metadata and sample data: to do this:

# set the working directory withthe downloaded metadata and sample data CSV files:

setwd("C:/PhD/Research/Data/Data_inventory/DOW")

# this is just some legacy code when I have more DOW files. It was to read them all in at once and putthem as dataframes in the envrionment (I could have just as easily done two calls of read.csv):

# read in csv file names in the working directory:

temp = list.files(pattern="*.csv") 

# then lapply read.csv function to the list of names:

dfs = lapply(temp, read.csv, fileEncoding="UTF-8-BOM")

# add the names of the file names as the names of the dfs in the list, minus the ".csv" (last 4 characters removed)

names(dfs) <- substr(temp,1,nchar(temp)-4) 

# map the items of the list to the global env.:

list2env(dfs,globalenv()) 

# rename dataframes to shorter names

sc<-Streams_Chemistry
sms<-Stream_Monitoring_Sites 

# filter chemistry for TP:

df.DOW_TP<-sc%>%filter(parameter_name=="PHOSPHORUS, TOTAL")

# I want to look at how many samples are there per unique site:

df.temp<-df.DOW_TP%>%group_by(site_id)%>%summarize(n=n())%>%
  arrange(desc(n))

# I want to see if these sites match up with a USGS gauge, but this is a pain in the ass!
# some sites might fall on the same stream that has gauge data, which would be easy to then use DA scaling to get the flows at the sample site
# but the workflow to do this is to do it by HAND (unless I want to automate it, which I could proablyfigure out but it wouldnt be perfect and would take time to learn how to do it)
# also, if the sample location is not on a stream, I would need to come up with a suitable nearby gauge, which again would be by HAND

# lets just start with making a df of DOW sites with number of TP samples > 20:

df.DOW_TP_sites<-left_join(df.temp,sms[,c(1,5,6)], by = 'site_id')%>%
  mutate(agency_cd = 'DOW Sample')%>%
  filter(n>20)

# and look at the map:

map.DOW.TP_sites<-df.DOW_TP_sites%>%
  drop_na(latitude,longitude)%>%
  st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)

mapview(map.DOW.TP_sites, zcol = 'agency_cd')

# lets add to this map the NWIS TP sites: to do this
# first combine the two dataframes:

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

####------------Step 4: CSI Import------------####

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

# set the WD:

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

df.CSI<-left_join(df.CSI, df.CSI_sites, by = 'Name') #%>% select(!Name)

# filter Sites for TP (also reformat the date column). Todo this:

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







# save.image(file = 'C:/PhD/Research/Code/CQ/Large_Scale_CQ_1.Rdata')

load('C:/PhD/Research/Code/CQ/Large_Scale_CQ_1.Rdata')








