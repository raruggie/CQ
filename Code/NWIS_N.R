# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

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


# TN pcode:

pcode<-'00600'
























#### NWIS Query ####

# this code adapted from NYS_site_ranking.R

# Download the metadata for sites with daily flow data. To do this:

# Use dataRetrieval::whatNWISdata to query USGS gauge sites with daily flow data in NY

df.NWIS.Q_sites<- whatNWISdata(stateCd = 'NY', parameterCd = "00060")

# write.csv(df.NWIS.Q_sites, "C:/PhD/CQ/Raw_Data/df.NWIS.Q_sites.csv", row.names=FALSE)

df.NWIS.Q_sites <- read.csv("Raw_Data/df.NWIS.Q_sites.csv", colClasses = c(site_no = "character"))

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

# From the sites with daily flow data, determine which sites also have over 99 discrete samples from one or more of the constituents of interests. To do this:

# use function I created (sourced) to get a dataframe of sites for just TN with #samples = 20 threshold:

df.NWIS.TN_sites<-fun.df.Pair_consit_flow(pcode, df.NWIS.Q_sites, n_samples = 20, state = 'NY')

# write.csv(df.NWIS.TN_sites, "Raw_Data/df.NWIS.TN_sites.csv", row.names=FALSE)

df.NWIS.TN_sites<-read.csv("Raw_Data/df.NWIS.TN_sites.csv", colClasses = c(site_no = "character"))

# download the raw daily flow data for these sites. To do this:
# readNWISdv can be used to download a single dataframe with all the raw flow data for a vector of gauge numbers (this takes a looong time):

df.NWIS.Q.for_TN<-readNWISdv(siteNumbers = df.NWIS.TN_sites$site_no, parameterCd = '00060', startDate = "", endDate = "", statCd = "00003")

# write.csv(df.NWIS.Q.for_TN, "Raw_Data/df.NWIS.Q.for_TN.csv", row.names=FALSE)

df.NWIS.Q.for_TN<-read.csv("Raw_Data/df.NWIS.Q.for_TN.csv", colClasses = c(site_no = "character"))

# download the raw discrete TN sample data:

df.NWIS.TN<-readNWISqw(siteNumbers = df.NWIS.TN_sites$site_no, parameterCd = pcode)
beepr::beep(8)
# write.csv(df.NWIS.TN, "Raw_Data/df.NWIS.TN.csv", row.names=FALSE)

df.NWIS.TN<-read.csv("Raw_Data/df.NWIS.TN.csv", colClasses = c(site_no = "character"))

#























#### Building CQ df ####

# join the TN with the flow df:

df.NWIS.TN_CQ<-left_join(df.NWIS.TN, df.NWIS.Q.for_TN, by=c("site_no"="site_no", "sample_dt"="Date"))

# remove observations where there are not CQ pairs:

df.NWIS.TN_CQ<-df.NWIS.TN_CQ%>%drop_na(X_00060_00003)

# take average of multiple samples on the same day at the same site:

df.NWIS.TN_CQ<-df.NWIS.TN_CQ%>%
  group_by(site_no, sample_dt)%>%
  summarise_at(vars(result_va, X_00060_00003), funs(mean(., na.rm=TRUE)))
  
# arrange by number of TN observations. To do this:

# first arrange the dataframe with the number of samples: 

temp<-df.NWIS.TN%>%group_by(site_no)%>%
  summarise(n=n())%>%
  arrange(desc(n))%>%
  mutate(n_sample_rank=rank(-n, ties.method='first'))

# merge this df with df.NWIS.TN_CQ and arrange by the new column:

df.NWIS.TN_CQ<-left_join(df.NWIS.TN_CQ,temp, by='site_no')%>%
  arrange(n_sample_rank)

# next step is to use Q yield to get better looking plot... idk if it will help but want to try also doesnt hurt to have the watershed areas as well. To do this:

# download the drainage areas from site metadata using readNWISsite

df.NWIS.TN_site_metadata<-readNWISsite(siteNumbers = unique(df.NWIS.TN_CQ$site_no))

# write.csv(df.NWIS.TN_site_metadata, "Raw_Data/df.NWIS.TN_site_metadata.csv", row.names=FALSE)

df.NWIS.TN_site_metadata<-read.csv("C:/PhD/CQ/Raw_Data/df.NWIS.TN_site_metadata.csv", colClasses = c(site_no = "character"))

# then select just the site number and DA column:

df.DA<-df.NWIS.TN_site_metadata%>%
  select(site_no, drain_area_va)

# finally merge with df.NWIS.TN_CQ and create a new Q column with area normalized flows (not worrying about units right now): 
# Note: will filter for NA in C and Q for breakpoint analysis in the next step as to keep the full list of sites with CQ pairs in this dataframe.
# Note: Some sites returned NA on draiange areas in readNWISsite, but I'll delinate anyways so I want the full list:

df.NWIS.TN_CQ<-left_join(df.NWIS.TN_CQ, df.DA, by = 'site_no')%>%
  mutate(Q_yield = X_00060_00003/drain_area_va)

#






















#### Fitting Breakpoints to CQ curves #### 

# create empty list to hold the results of the segmented function (which does the breakpoint analysis)

l_Seg<-list()

# create a matrix to hold the results of the davies test, which determines if a two slope model is warrented over a single slope model:

davies.test.matrix<-NULL

# create a new dataframe of only paired CQ observations (such that the breakpoit analysis function runs smoothly) (I didnt want to do this in loop for some reason, I cant remeber why but it wouldnt work):

df.NWIS.TN_CQ_for_BP<-df.NWIS.TN_CQ%>%
  drop_na(result_va, Q_yield)

# create a vector of ordered unique site names:

temp.loop<-sort(unique(df.NWIS.TN_CQ_for_BP$site_no))

# test i for for loop building:

# i<-4

# loop through the sites:

for (i in seq_along(temp.loop)){
  
  tryCatch({
    
    # print the site name for loop debugging:
    
    print(i)
    print(temp.loop[i])
    
    # create a dataframe that will work with segmented. To do this: 
      # filter for the site the loop is in
      # add log transformed C and Q columns, as well as duplicated columns for renamed C and Q
      # filter for real log C and Q values so breakpoint analysis works smoothly:
    
    df<-df.NWIS.TN_CQ_for_BP%>%
      filter(site_no == temp.loop[i])%>%
      mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
      filter(is.finite(log_C))%>%
      filter(is.finite(log_Q))
    
    # build a single slope lm for log C and Q. Tis model is also used inthe breakpoint analysis inthenext step:
    
    m<-lm(log_C~log_Q, df)
    
    # perform breakpoint regression:
    
    m_seg<-segmented(obj = m, npsi = 1)
    
    # perform davies test for constant linear predictor:
    # the results are saved as a string with the the site name and true/false:
    
    x<-paste(temp.loop[i], '-', davies.test(m)$p.val<0.05)
    
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
      result_df<-fit%>%mutate(site = temp.loop[i], Date = df$sample_dt, n = df$n, Q_real = df$Q, C = df$C, I1 = inter$Est.[1], I2 = inter$Est.[2], Slope1 = s[1,1], Slope2 = s[2,1], BP = bp)
    } else{
      result_df<-fit%>%mutate(site = temp.loop[i], Date = df$sample_dt, n = df$n, Q_real = df$Q, C = df$C, I1 = NA, I2 = NA, Slope1 = NA, Slope2 = NA, BP = NA)
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

# now combine the list of dfs of the breakpoint analysis results (with fitted values, intercepts and slopes) into a single df:

df_Seg<-bind_rows(l_Seg)

# merge this dataframe with the davies test result dfs to add the BP_yes column, use replace and the BP yes column with a conditional statement to set the breakpoint Q and C column rows to NA, as to not plot the segmeneted line if davies test was false:

df_Seg<-df_Seg%>%
  left_join(., df.davies, by = 'site')%>%
  mutate(across(c(1,2), ~replace(., BP_yes == 'FALSE', NA)))

# add a column for number of samples:

df_Seg<-left_join(df_Seg, temp[,c(1,3)], by = c('site'='site_no'))%>%
  arrange(n_sample_rank)

# ready to plot:

p<-ggplot(df_Seg, aes(x = log(Q_real), y = log(C)))+
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'tomato')+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# p


# one last thing: lets look at a map of these:

map.NWIS.TN_sites<-df.NWIS.TN_site_metadata%>%
  rename(longitude=8,latitude=7)%>%
  drop_na(latitude,longitude)%>%
  st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)%>%
  mutate(site_id= paste(site_no, station_nm), n = NA, .before = 1)%>%
  select(c(1:3))

mapview(map.NWIS.TN_sites)



























#### Tradeoff matrix ####

# build a matrix of the number of TN samples as a function of watershed size. To do this:

# I already have a df with paired CQ observations and DA, just need to get the distinct sites:

TN_sites<-df.NWIS.TN_CQ%>%distinct(site_no, .keep_all = T)

# set up df to populate:

m<-data.frame(Min_num_samples  = c(20,50,75,100,200), '25' = NA, '50' = NA, '100' = NA, '150'=NA, '250'=NA, '500'=NA, '1000'=NA, 'Unlimited'=NA)

# set up variables for thresholds for number of samples and DA size:

n_sam<- c(20,50,75,100,200)-1

min_DA<- c(25,50,100,150,250,500,1000)

# loop through number of samples (rows):

# i<-1

for (i in seq_along(n_sam)){
  
  temp.i<-TN_sites%>%filter(n>n_sam[i])
  
  m$Unlimited[i]<-dim(temp.i)[1]

  # j<-2
  
  # loop through the size of the DA (columns)
  for(j in  seq_along(min_DA)){
    
    temp.j<-temp.i%>%filter(drain_area_va<=min_DA[j])
    
    m[i,j+1]<-dim(temp.j)[1]
    
  }
  
}

m

#



















####~~~~ Test to see if N sites are in Gauges 2 (Not run) ~~~~####

# read in the gauges 2 database (I am forgoing the datalayers workflow for now since I may have a great set of predictors in this database): to do this:

# read in all sheets using function

l.G2 <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm.xlsx")

# remove the last element (does notcomtain useful info)

l.G2[[27]]<-NULL

# and convert to a single df

df.G2<-reduce(l.G2, full_join, by = "STAID")

# filter the gauges2 to the NYS TN CQ sites:

df.G2<-df.G2%>%filter(STAID %in% TN_sites$site_no)

# only 68 of the orginal 102 TN sites are in gauges 2

























####~~~~ Filter down to final site list for manuscript ~~~~####

# filter the raw CQ df based on date and to get the desired sites:
# note: this was initnal employed because of the filter for sites withsuccessful delineations
# nut I am not including that requirment any more for the CQ manuscript
# but I am still keeping the 2001 requirment because gauges 2 has NLCD 2006, which I amsati

temp<-df.NWIS.TN_CQ%>%filter(year(sample_dt) >= 2001)

unique(temp$site_no)

# 39 sites

# map of this new set of sites:

df.NWIS.TN_site_metadata%>%
  filter(site_no %in% temp$site_no)%>%
  st_as_sf(.,coords=c('dec_long_va','dec_lat_va'), crs = 4326, remove = FALSE)%>%
  mapview(.)

# remove sites below a certain latitude to get rid of long island:

temp1<-filter(df.NWIS.TN_site_metadata, dec_lat_va >40.9364)

# look at map:

temp1%>%
  st_as_sf(.,coords=c('dec_long_va','dec_lat_va'), crs = 4326, remove = FALSE)%>%
  mapview(.)

# now filter final df:

temp<-temp%>%filter(site_no %in% temp1$site_no)

unique(temp$site_no)

# still 39 sites

# filter on gauges 2:

temp<-filter(temp, site_no %in% df.G2$STAID)

unique(temp$site_no)

# 23 sites

# map:

df.NWIS.TN_site_metadata%>%
  filter(site_no %in% temp$site_no)%>%
  st_as_sf(.,coords=c('dec_long_va','dec_lat_va'), crs = 4326, remove = FALSE)%>%
  mapview(.)

# number of sites ingauges 2 and not on LI (i.e. without post 2001 filter):

length(unique((temp1%>%filter(site_no%in% df.G2$STAID))$site_no))

#



# save.image(file = 'Processed_Data/NWIS_N.Rdata')
load('Processed_Data/NWIS_N.Rdata')

#



