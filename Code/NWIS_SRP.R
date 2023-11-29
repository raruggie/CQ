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

# use function I created (sourced) to get a dataframe of sites for just TP with #samples = 20 threshold:

df.NWIS.SRP_sites<-fun.df.Pair_consit_flow('00660', df.NWIS.Q_sites, n_samples = 20, state = 'NY')
beepr::beep(8)
write.csv(df.NWIS.SRP_sites, "Raw_Data/df.NWIS.SRP_sites.csv", row.names=FALSE)

df.NWIS.SRP_sites<-read.csv("Raw_Data/df.NWIS.SRP_sites.csv", colClasses = c(site_no = "character"))

# download the raw daily flow data for these sites. To do this:
# readNWISdv can be used to download a single dataframe with all the raw flow data for a vector of gauge numbers (this takes a looong time):

df.NWIS.Q.for_SRP<-readNWISdv(siteNumbers = df.NWIS.SRP_sites$site_no, parameterCd = '00060', startDate = "", endDate = "", statCd = "00003")

# write.csv(df.NWIS.Q.for_SRP, "Raw_Data/df.NWIS.Q.for_SRP.csv", row.names=FALSE)

df.NWIS.Q.for_SRP<-read.csv("Raw_Data/df.NWIS.Q.for_SRP.csv", colClasses = c(site_no = "character"))

# download the raw discrete TP sample data:

df.NWIS.SRP<-readNWISqw(siteNumbers = df.NWIS.SRP_sites$site_no, parameterCd = '00600')

# write.csv(df.NWIS.SRP, "Raw_Data/df.NWIS.SRP.csv", row.names=FALSE)

df.NWIS.SRP<-read.csv("Raw_Data/df.NWIS.SRP.csv", colClasses = c(site_no = "character"))

#























#### Building CQ df ####

# join the TP with the flow df:

df.NWIS.SRP_CQ<-left_join(df.NWIS.SRP, df.NWIS.Q.for_SRP, by=c("site_no"="site_no", "sample_dt"="Date"))

# remove observations where there are not CQ pairs:

df.NWIS.SRP_CQ<-df.NWIS.SRP_CQ%>%drop_na(X_00060_00003)

# take average of multiple samples on the same day at the same site:

df.NWIS.SRP_CQ<-df.NWIS.SRP_CQ%>%
  group_by(site_no, sample_dt)%>%
  summarise_at(vars(result_va, X_00060_00003), funs(mean(., na.rm=TRUE)))
  
# arrange by number of TP observations. To do this:

# first arrange the dataframe with the number of samples: 

temp<-df.NWIS.SRP%>%group_by(site_no)%>%
  summarise(n=n())%>%
  arrange(desc(n))%>%
  mutate(n_sample_rank=rank(-n, ties.method='first'))

# merge this df with df.NWIS.TP_CQ and arrange by the new column:

df.NWIS.SRP_CQ<-left_join(df.NWIS.SRP_CQ,temp, by='site_no')%>%
  arrange(n_sample_rank)

# next step is to use Q yield to get better looking plot... idk if it will help but want to try also doesnt hurt to have the watershed areas as well. To do this:

# download the drainage areas from site metadata using readNWISsite

df.NWIS.SRP_site_metadata<-readNWISsite(siteNumbers = unique(df.NWIS.SRP_CQ$site_no))

# then select just the site number and DA column:

df.DA<-df.NWIS.SRP_site_metadata%>%
  select(site_no, drain_area_va)

# finally merge with df.NWIS.SRP_CQ and create a new Q column with area normalized flows (not worrying about units right now): 
# Note: will filter for NA in C and Q for breakpoint analysis in the next step as to keep the full list of sites with CQ pairs in this dataframe.
# Note: Some sites returned NA on draiange areas in readNWISsite, but I'll delinate anyways so I want the full list:

df.NWIS.SRP_CQ<-left_join(df.NWIS.SRP_CQ, df.DA, by = 'site_no')%>%
  mutate(Q_yield = X_00060_00003/drain_area_va)

#






















#### Fitting Breakpoints to CQ curves #### 

# create empty list to hold the results of the segmented function (which does the breakpoint analysis)

l_Seg<-list()

# create a matrix to hold the results of the davies test, which determines if a two slope model is warrented over a single slope model:

davies.test.matrix<-NULL

# create a new dataframe of only paired CQ observations (such that the breakpoit analysis function runs smoothly) (I didnt want to do this in loop for some reason, I cant remeber why but it wouldnt work):

df.NWIS.SRP_CQ_for_BP<-df.NWIS.SRP_CQ%>%
  drop_na(result_va, Q_yield)

# create a vector of ordered unique site names:

temp.loop<-sort(unique(df.NWIS.SRP_CQ_for_BP$site_no))

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
    
    df<-df.NWIS.SRP_CQ_for_BP%>%
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

p


# one last thing: lets look at a map of these:

map.NWIS.SRP_sites<-df.NWIS.SRP_site_metadata%>%
  rename(longitude=8,latitude=7)%>%
  drop_na(latitude,longitude)%>%
  st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)%>%
  mutate(site_id= paste(site_no, station_nm), n = NA, .before = 1)%>%
  select(c(1:3))

mapview(map.NWIS.SRP_sites)



























#### Tradeoff matrix ####

# build a matrix of the number of TP samples as a function of watershed size. To do this:

# I already have a df with paired CQ observations and DA, just need to get the distinct sites:

SRP_sites<-df.NWIS.SRP_CQ%>%distinct(site_no, .keep_all = T)

# set up df to populate:

m<-data.frame(Min_num_samples  = c(20,50,75,100,200), '25' = NA, '50' = NA, '100' = NA, '150'=NA, '250'=NA, '500'=NA, '1000'=NA, 'Unlimited'=NA)

# set up variables for thresholds for number of samples and DA size:

n_sam<- c(20,50,75,100,200)-1

min_DA<- c(25,50,100,150,250,500,1000)

# loop through number of samples (rows):

# i<-1

for (i in seq_along(n_sam)){
  
  temp.i<-SRP_sites%>%filter(n>n_sam[i])
  
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

# filter the gauges2 to the NYS TP CQ sites:

df.G2<-df.G2%>%filter(STAID %in% SRP_sites$site_no)

# only 48 of the orginal 61  sites are in gauges 2





















# save.image(file = 'Processed_Data/NWIS_SRP.Rdata')
load('Processed_Data/NWIS_SRP.Rdata')


#




















#### Delinating Watersheds (NOT RUN) ####

# # delinate the 137 NWIS sites:
# 
# l.SS_WS.NWIS<-lapply(seq_along(df.NWIS.TN_site_metadata$site_no), \(i) Delineate(df.NWIS.TP_site_metadata$dec_long_va[i], df.NWIS.TP_site_metadata$dec_lat_va[i]))
# 
# names(l.SS_WS.NWIS)<-df.NWIS.TP_site_metadata$site_no
# 
# # save(l.SS_WS.NWIS, file = 'C:/PhD/CQ/Downloaded_Data/l.SS_WS.NWIS.Rdata')
# 
# df.sf.NWIS<-fun.l.SS_WS.to.sfdf(l.SS_WS.NWIS)
# 
# # save(df.sf.NWIS, file = 'C:/PhD/CQ/Processed_Data/df.sf.NWIS.Rdata')
# 
# load('C:/PhD/CQ/Processed_Data/df.sf.NWIS.Rdata')
# 
# # calculate DA using vect:
# 
# vect.NWIS<-vect(df.sf.NWIS)
# 
# vect.NWIS$area_KM2<-expanse(vect.NWIS, unit="km")
# 
# # add DA in mi2 to df.sf:
# 
# df.sf.NWIS$area_sqmi<-expanse(vect.NWIS, unit="km")*0.386102
# 
# # add NWIS tabulated area to df sf: first download the metadata for the sites using readNWISsite, then merge drain_area_va column to sfdf:
# 
# temp<-readNWISsite(df.sf.NWIS$Name)
# 
# df.sf.NWIS<-left_join(df.sf.NWIS, temp[,c(2,30)], by = c('Name'='site_no'))
# 
# # calcuate the percent error of the delination:
# 
# df.sf.NWIS$Delination_Error<-(df.sf.NWIS$area_sqmi-df.sf.NWIS$drain_area_va)/df.sf.NWIS$drain_area_va
# 
# # I want to look at the delimaitons that are on the cusp of not working:
# 
# temp<-df.sf.NWIS[,c(1,107)]%>%
#   filter(abs(Delination_Error)>.02)
# 
# # mapview(temp)
# 
# # after playing with the numbers, 2% seems like a good threshold for non-delineations:
# # subset df.sf to those sites that are considered to delineate correctly:
# 
# df.sf.NWIS.keep<-df.sf.NWIS%>%filter(abs(Delination_Error)<=.02)
# 
# # note I suspect that snapping the lat long to nhd would improve the delineation success notgoing to do that now
# 





















####~~~~ Filter down to final site list for manuscript ~~~~####

# filter the raw CQ df based on date and to get the desired sites:
# note: this was initnal employed because of the filter for sites withsuccessful delineations
# nut I am not including that requirment any more for the CQ manuscript
# but I am still keeping the 2001 requirment because gauges 2 has NLCD 2006, which I amsati

temp<-df.NWIS.TN_CQ%>%filter(year(sample_dt) >= 2001)

# 39 sites

# map of this new set of sites:

df.NWIS.TN_site_metadata%>%
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

length(unique((temp1%>%
                 filter(site_no%in% df.G2$STAID))$site_no))

#






####~~~~ Filter down to final site list for manuscript ~~~~####

# filter the raw CQ df based on date and to get the desired sites:
# note: this was initnal employed because of the filter for sites withsuccessful delineations
# nut I am not including that requirment any more for the CQ manuscript
# but I am still keeping the 2001 requirment because gauges 2 has NLCD 2006, which I amsati

temp<-df.NWIS.SRP_CQ%>%filter(year(sample_dt) >= 2001)

# 39 sites

# map of this new set of sites:

df.NWIS.SRP_site_metadata%>%
  st_as_sf(.,coords=c('dec_long_va','dec_lat_va'), crs = 4326, remove = FALSE)%>%
  mapview(.)

# remove sites below a certain latitude to get rid of long island:

temp1<-filter(df.NWIS.SRP_site_metadata, dec_lat_va >40.9364)

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

# 9 sites

# map:

df.NWIS.SRP_site_metadata%>%
  filter(site_no %in% temp$site_no)%>%
  st_as_sf(.,coords=c('dec_long_va','dec_lat_va'), crs = 4326, remove = FALSE)%>%
  mapview(.)

#























############################################################
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
####~~~~~~~~~~~~~~~~~Watershed Attributes~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~(for correlation, eventually)~~~~~~~~~~~####
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
############################################################


























####~~~~ Gauges 2 (Not run) ~~~~####

# read in the gauges 2 database (I am forgoing the datalayers workflow for now since I may have a great set of predictors in this database): to do this:

# read in all sheets using function

l.G2 <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm.xlsx")

# remove the last element (does notcomtain useful info)

l.G2[[27]]<-NULL

# and convert to a single df

df.G2<-reduce(l.G2, full_join, by = "STAID")

# filter the gauges2 to the NYS TP CQ sites:

df.G2<-df.G2%>%filter(STAID %in% TP_sites$site_no)

# only 89 of the orginal 137 sites are in gauges 2

# lets see what the tradeoff matrix looks like if only these 89 are included:

TP_sites.2<-df.NWIS.TP_CQ%>%
  filter(site_no %in% df.G2$STAID)%>%
  distinct(site_no, .keep_all = T)

m<-data.frame(Min_num_samples  = c(20,50,75,100,200), '25' = NA, '50' = NA, '100' = NA, '150'=NA, '250'=NA, '500'=NA, '1000'=NA, 'Unlimited'=NA)

n_sam<- c(20,50,75,100,200)-1

min_DA<- c(25,50,100,150,250,500,1000)

# i<-1

for (i in seq_along(n_sam)){
  
  temp.i<-TP_sites.2%>%filter(n>n_sam[i])
  
  m$Unlimited[i]<-dim(temp.i)[1]
  
  # j<-2
  
  for(j in  seq_along(min_DA)){
    
    temp.j<-temp.i%>%filter(drain_area_va<=min_DA[j])
    
    m[i,j+1]<-dim(temp.j)[1]
    
  }
  
}

m

# We didnt end up using this tradeoff matrix
# the issue is that the attributs that comewiththe streamstatsdownloads I dont think are accurate
# (I thinkthere are lots of extra zeros)
# this makes me want to use the gauges 2 sites for the correlaitons/regressions with CQ parameters

# I am writing post going through the datalayer downloads and processing below


























####~~~~ Data Layers ~~~~####

#### CDL: ####

# note when I ran ths CDL code block I used the 103 sites thatwere apartofthe set, this is prior to filtering based onthe CDLdate... so I am just going to filter the resulting CDL df at the end. Just note that what is saved is of the 103 sites, but that the code doesnt reflect this, i.e. I just ran the filter anddidnt keep the code because df.NWIS.TP.keep would be down to 56 if I did it right 

# download:

rast.NWIS.CDL.2020 <- GetCDLData(aoi = df.sf.NWIS, year = "2020", type = "b", tol_time = 1000)
rast.NWIS.CDL.2008 <- GetCDLData(aoi = df.sf.NWIS, year = "2008", type = "b", tol_time = 1000)

# convert to rast:

rast.NWIS.CDL.2020<-rast(rast.NWIS.CDL.2020)
rast.NWIS.CDL.2008<-rast(rast.NWIS.CDL.2008)

# check CRS of both years:

crs(rast.NWIS.CDL.2020)
crs(rast.NWIS.CDL.2008)

# plot

# plot(rast.NWIS.CDL.2020)
# plot(rast.NWIS.CDL.2008)

# reproject to sample watershed vector data to match raster data:

vect.NWIS<-vect(df.sf.NWIS) # need to first ininalize vect since sometimes reading in rdata file (terra issue)

vect.NWIS.proj<-terra::project(vect.NWIS, crs(rast.NWIS.CDL.2008))

# extract frequency tables for each sample watershed

l.NWIS.CDL <- terra::extract(rast.NWIS.CDL.2020, vect.NWIS.proj, table,ID=FALSE)[[1]]
l.NWIS.CDL.2008 <- terra::extract(rast.NWIS.CDL.2008, vect.NWIS.proj, table,ID=FALSE)[[1]]

# save:

# save(l.NWIS.CDL, file='Processed_Data/l.NWIS.CDL.Rdata')
# save(l.NWIS.CDL.2008, file='Processed_Data/l.NWIS.CDL.2008.Rdata')

load('Processed_Data/l.NWIS.CDL.Rdata')
load('Processed_Data/l.NWIS.CDL.2008.Rdata')

# convert resulting list of tables to list of dfs

l.NWIS.CDL<-lapply(l.NWIS.CDL, as.data.frame)
l.NWIS.CDL.2008<-lapply(l.NWIS.CDL.2008, as.data.frame)

# the next step is to join the CDL key with these dfs, but first:

# aggregate CDL: to do this:
# looking at the CDL legend (linkdata), I determined the following:

Ag<-c(1:6,10:14,21:39,41:61,66:72,74:77,204:214,216:227,229:250,254)
Pasture<-c(176)
Forest<-c(63,141:143)
Developed<-c(82,121:124)
Water<-c(83,111)
Wetlands_all<-c(87,190,195)
Other<-c(64,65,88,112,131,152) # Shrub<-c(64,152) Barren<-c()

l <- tibble::lst(Ag,Pasture,Forest,Developed,Water,Wetlands_all,Other)

reclass_CDL<-data.frame(lapply(l, `length<-`, max(lengths(l))))%>%
  pivot_longer(cols = everything(), values_to = 'MasterCat',names_to = 'Crop')%>%
  drop_na(MasterCat)

# left join each df in the list to the CDL legend key, as well as calcuate the pland:
  
l.NWIS.CDL<-lapply(l.NWIS.CDL, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., reclass_CDL, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes
l.NWIS.CDL.2008<-lapply(l.NWIS.CDL.2008, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., reclass_CDL, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes

# set names of list:

names(l.NWIS.CDL)<-df.sf.NWIS$Name # note doesnt work with the workflow set up to filter on 2008 limiter
names(l.NWIS.CDL.2008)<-df.sf.NWIS$Name # note doesnt work with the workflow set up to filter on 2008 limiter

# combine list into single df:

df.NWIS.CDL<-bind_rows(l.NWIS.CDL, .id = 'Name')
df.NWIS.CDL.2008<-bind_rows(l.NWIS.CDL.2008, .id = 'Name')

# remove potential for one of the MasterCats in the CDL to be empty, which is messing with the pivot_wider below:

df.NWIS.CDL<-filter(df.NWIS.CDL, Crop != '')
df.NWIS.CDL.2008<-filter(df.NWIS.CDL.2008, Crop != '')

# pivot wider:

# if using the CDL linkdata, use this:

# df.NWIS.CDL<- pivot_wider(df.NWIS.CDL[,-2], names_from = Crop, values_from = Freq)
# df.NWIS.CDL.2008<- pivot_wider(df.NWIS.CDL.2008[,-2], names_from = Crop, values_from = Freq)

# if using the reclass_CDL data, use this:

df.NWIS.CDL<-df.NWIS.CDL[,-2]%>%
  group_by(Name, Crop)%>%
  summarise(Freq=sum(Freq, na.rm = T))%>%
  pivot_wider(., names_from = Crop, values_from = Freq)

df.NWIS.CDL.2008<-df.NWIS.CDL.2008[,-2]%>%
  group_by(Name, Crop)%>%
  summarise(Freq=sum(Freq, na.rm = T))%>%
  pivot_wider(., names_from = Crop, values_from = Freq)

# check to see if add up to 100%:

sort(rowSums(df.NWIS.CDL[,-1], na.rm = T))
sort(rowSums(df.NWIS.CDL.2008[,-1], na.rm = T))

# one of the CDL 2020 is low:

which.min(rowSums(df.NWIS.CDL[,-1], na.rm = T))

# this is a long island site, so doesnt matter

# other than that site looks good. 

# Done with CDL

# actually I want to come up with a workflow to get all years of CDL for the watersheds,
# butthat isgoing totake a long timeto run...

# but first I'm going to look at the differences in the CDL between 2008 and 2020


























#### NLCD ####

# the CDL and NLCD differ in that the CDL combined Pasture and grassland while the nLCD does not
# thus, even though the CDL can be agrgated to et close to the NLCD, I want to also have the nLCD
# to the difference.

# note: the 53 sites post filtering based onthe CDLdate are used in this workflow

# I am going to run this workflow for NLCD 2019 and NLCD 2001 to see if any sites had major changes in land use:

# download: to do this:
# NLCD is not working when trying to download based on the entire polygon df, so going to use lapply to download individually:

l.rast.NWIS.NLCD.2019 <- lapply(seq_along(df.sf.NWIS.keep$Name), \(i) get_nlcd(template = st_cast(df.sf.NWIS.keep, "MULTIPOLYGON")[i,], label = as.character(i), year = 2019))
l.rast.NWIS.NLCD.2001 <- lapply(seq_along(df.sf.NWIS.keep$Name), \(i) get_nlcd(template = st_cast(df.sf.NWIS.keep, "MULTIPOLYGON")[i,], label = as.character(i), year = 2001))

# save(l.rast.NWIS.NLCD.2019, file='Downloaded_Data/l.rast.NWIS.NLCD.2019.Rdata')
# save(l.rast.NWIS.NLCD.2001, file='Downloaded_Data/l.rast.NWIS.NLCD.2001.Rdata')

# convert to SpatRasters:

l.rast.NWIS.NLCD.2019<-lapply(l.rast.NWIS.NLCD.2019, rast)
l.rast.NWIS.NLCD.2001<-lapply(l.rast.NWIS.NLCD.2001, rast)

# plot

# plot(rast.NWIS.NLCD.2019)

# see if 2001 and 2019 crs are the same:

crs(l.rast.NWIS.NLCD.2019[[2]])
crs(l.rast.NWIS.NLCD.2001[[2]])

# reproject to sample watershed vector data to match raster data:

vect.NWIS<-vect(df.sf.NWIS.keep) # need to first ininalize vect since sometimes reading in rdata file (terra issue)

vect.NWIS.proj<-terra::project(vect.NWIS, crs(l.rast.NWIS.NLCD.2019[[1]]))

# extract frequency tables for each sample watershed

system.time({l.NWIS.NLCD.2019 <- lapply(seq_along(l.rast.NWIS.NLCD.2019), \(i) terra::extract(l.rast.NWIS.NLCD.2019[[i]], vect.NWIS.proj[i], ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2)))})
system.time({l.NWIS.NLCD.2001 <- lapply(seq_along(l.rast.NWIS.NLCD.2001), \(i) terra::extract(l.rast.NWIS.NLCD.2001[[i]], vect.NWIS.proj[i], ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2)))})

# save(l.NWIS.NLCD.2019, file = 'Processed_Data/l.NWIS.NLCD.2019.Rdata')
# save(l.NWIS.NLCD.2001, file = 'Processed_Data/l.NWIS.NLCD.2001.Rdata')

load("Processed_Data/l.NWIS.NLCD.2019.Rdata")
load("Processed_Data/l.NWIS.NLCD.2001.Rdata")

# reclassify: to do this:

# adjust the NLCD reclassify legend to match the CDL reclassify df made above:

legend.NWIS<-legend; legend.NWIS$Class3[c(13,17)]<-'Pasture';legend.NWIS$Class3[14]<-'Other';legend.NWIS$Class3[1]<-'Water';legend.NWIS$Class3[c(19,20)]<-'Wetlands_all'

sort(unique(legend.NWIS$Class3))==sort(unique(reclass_CDL$Crop))

# looks good

# reclassify the NLCD using this new legend and clean up the dataframe from the next step

l.NWIS.NLCD.2019<-lapply(l.NWIS.NLCD.2019, \(i) left_join(as.data.frame(i), legend.NWIS%>%select(Class, Class3), by = 'Class')%>%mutate(Class = Class3)%>%select(-Class3))
l.NWIS.NLCD.2001<-lapply(l.NWIS.NLCD.2001, \(i) left_join(as.data.frame(i), legend.NWIS%>%select(Class, Class3), by = 'Class')%>%mutate(Class = Class3)%>%select(-Class3))

# pivot_wider the df in the lists and add a Name column for the site:

l.NWIS.NLCD.2019<-lapply(seq_along(l.NWIS.NLCD.2019), \(i) l.NWIS.NLCD.2019[[i]]%>%group_by(Class)%>%summarise(Freq = sum(Freq))%>%pivot_wider(names_from = Class, values_from = Freq)%>%mutate(Name = df.sf.NWIS.keep$Name[i], .before = 1)%>%as.data.frame(.))
l.NWIS.NLCD.2001<-lapply(seq_along(l.NWIS.NLCD.2001), \(i) l.NWIS.NLCD.2001[[i]]%>%group_by(Class)%>%summarise(Freq = sum(Freq))%>%pivot_wider(names_from = Class, values_from = Freq)%>%mutate(Name = df.sf.NWIS.keep$Name[i], .before = 1)%>%as.data.frame(.))

# bind the lists into a single dataframe
# note some of the sites have different length dtaframes because they didnt have all the same number of NLCD classes. When binding rows this will give a dataframe of the maximum length and put NAs for sites where there wasn't a column: 

df.NWIS.NLCD.2019<-bind_rows(l.NWIS.NLCD.2019)
df.NWIS.NLCD.2001<-bind_rows(l.NWIS.NLCD.2001)

# create a list of these two dataframes:

l.NWIS.NLCD<-list(df.NWIS.NLCD.2019,df.NWIS.NLCD.2001)%>%purrr::set_names(c('2019','2001'))

# save thisprocessed list:

# save(l.NWIS.NLCD, file='Processed_Data/l.NWIS.NLCD.Rdata')

load('Processed_Data/l.NWIS.NLCD.Rdata')

#

























#### NED ####

# download:

DEM.NWIS<-get_ned(df.sf.NWIS, label = '2') # already SpatRaster!

writeRaster(DEM.NWIS, file='Downloaded_Data/DEM.NWIS.tif', overwrite=TRUE)

# plot(DEM.NWIS)

# extract elevation metrics over each sample watershed: to do this:
# build a function with multiple functions:

f <- function(x, na.rm = T) {
  c(mean=mean(x, na.rm = na.rm),
    range=max(x, na.rm = na.rm)-min(x, na.rm = na.rm),
    sd=sd(x, na.rm = na.rm)
  )
}

# reproject NWIS basins to DEM crs:

vect.NWIS<-vect(df.sf.NWIS)

vect.NWIS.proj<-terra::project(vect.NWIS, crs(DEM.NWIS))

# extract the metrics over each watershed using the function above:

df.NWIS.DEM <- as.data.frame(terra::extract(DEM.NWIS, vect.NWIS.proj, f))

# set the names of the df:

names(df.NWIS.DEM)<-c('Name', 'Elev_Avg', 'Elev_Range', 'Elev_SD')

# set the names of the sites:

df.NWIS.DEM$Name<-df.sf.NWIS$Name

# finally save the df:

# save(df.NWIS.DEM, file= 'Processed_Data/df.NWIS.DEM.Rdata')

load('Processed_Data/df.NWIS.DEM.Rdata')

# done with DEM

























#### Climate (not run) ####

# steps:

# Loop through the three metrics and metric specific operations.
# for each metric, perform:

# 1) download all gridded climate data for a metric (e.g.: precip)
# 2) extract daily mean of the metric for each subbasin into df
# 3) reformat (pivot longer),convert the date, and calcualte annual statistics based on the metric (e.g. precip is annual totals)

# the definition of:
# pr = total annual precio
# tmmn = minimum of all annual daily minimum temperature 
# tmmx = maximum of all annual daily maximum temperature

vars<-c('pr', 'tmmn', 'tmmx')

vars_funs<-c(`sum`, `min`, `max`)

# set up a df for appending into each loop:

df.NWIS.Climate<-data.frame(ID = 1:dim(df.NWIS.CDL)[1])

i<-1

# loop through the climate variables in vars:

for (i in seq_along(vars)){
  
  # download a Spatraster for the entire bounding box of all the watersheds of the climate variable data for the loops iteration
  # each layer of the spatraster is a day for the climate variable in thedaterange provided inthe function call:
  
  climate<-getGridMET(vect.NWIS, varname = vars[i], startDate = "2016-01-01", endDate = "2019-12-31")[[1]]
  
  # reproject the watersheds to the crs of the climate variables raster:
  
  vect.NWIS.proj<-terra::project(vect.NWIS, crs(climate))
  
  # extract the daily mean of the climate variables data over the watersheds:
  # thisreturns a dataframe thefirstcolumn is the site (justnumbered 1-56 right now) and the other columns are the daily values of the climate vairable (taking the mean ofthe raster cells over the watershed): 
  
  climate <- terra::extract(climate, vect.NWIS.proj, mean)
  
  # reformat, convert dates, and calcualte annual stats:
  
  # pivot Longer: each of the columns for the daily values of the climate variable starts withthevariable name and an underscore, so using that pattern to make a long dataframe 
  # format the dates: the date is contained withinthe column that was priviously pivoted into a single column. so reformating this column and turing into a date:
  # grouping by the site and year, summarize the daily values into an annualvalue based onthe operation for the climatevariable. for pricep, sum is used, for min temp min is used, and formax temp max is used!:
  # rename the reulting year(Date) column to a nicer name:
  # renaming the values ofthe year column: add the climate variable to the year,this pairs withthe next andfinal step:
  # pivot wider to get unique columns for the years of the cimlatevariable (here its 2016-2019):
  
  climate1<-climate%>%
    pivot_longer(cols= starts_with(paste0(vars[i],"_")), names_to = 'Date',values_to = "Value")%>%
    mutate(Date=as.Date(str_replace(Date, paste0(vars[i],"_"), "")))%>%
    group_by(year(Date),ID)%>%
    summarize(Annual = vars_funs[[i]](Value))%>%
    rename(Year=1)%>%
    mutate(Year=paste0(vars[i],Year))%>%
    pivot_wider(names_from = Year, values_from = Annual)
  
  # merge the climate variables data to a common dataframe (initalized prior to loop).
  # New columns are added for each site:
  
  df.NWIS.Climate<-left_join(df.NWIS.Climate, climate1, by = 'ID')
  
}

# rename the sites:

df.NWIS.Climate[,1]<-df.NWIS.DEM$Name

# rename the sites:

names(df.NWIS.Climate)[1]<-'Name'

# save df:

save(df.NWIS.Climate, file = 'Processed_Data/df.NWIS.Climate.Rdata')

#























##### Determine if CAFO is in watershed ####

# read in CAFO locations and make sf dataframe

CAFOs<-read.csv('Raw_data/Copy of CAFO-list-from-NYSDEC March 2021 with lat long values.csv')%>%
  drop_na(longitude)%>%
  drop_na(latitude)%>%
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326)

# add CAFOS to watershed map:

# mapview(df.sf.NWIS.keep)+
#   mapview(CAFOs, zcol = 'SIZE')

# this file maps differently than the map on the DE website, but just going
# to proceed with this file for now

# determine how many CAFOs intersect each watershed:

df.sf.NWIS.keep$CAFO_count <- lengths(st_intersects(df.sf.NWIS.keep, CAFOs))

# make a map colored by CAFO count:

# mapview(df.sf.NWIS.keep, zcol = 'CAFO_count')

#























#####################################################
####~~~~~~ Determine a final set of sites ~~~~~~~####
#####################################################






















####~~~~ Data Layers ~~~~####

# filter the raw CQ df based on date and to get the desired sites:

temp<-df.NWIS.TP_CQ%>%filter(year(sample_dt) >= 2001)

# then filter the sites:

df.datalayers<-df.sf.NWIS.keep%>%filter(Name %in% unique(temp$site_no))

# map of this new set of sites:

# mapview(df.datalayers)

# and remove the Long island sites:

df.datalayers<-df.datalayers%>%filter(!Name %in% c("01304000", "01305000", "01304500"))

#



















###################################################
 #### ~~~~ Now filtering the data layers ~~~~ ####
 #### ~~~~ sites to those in gauges 2    ~~~~ ####
###################################################






















####~~~~ Gauges 2 ~~~~####

# how many of the data layers sites are in gauges 2: to do this:

# load in gauges 2 sites (predictor set edited_again):
# (may need to redo source the functions code):

l.G2 <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm_EDITED_again.xlsx")

# reduce the sheets (list elements) into a single df:

df.G2<-reduce(l.G2, full_join, by = "STAID")

# select the first column andthe numeric ones after:

df.G2<-select(df.G2, STAID | where(is.numeric))

# filter the gauges2 to the sites IDed above:

df.G2<-df.G2%>%filter(STAID %in% df.datalayers$Name)
  
# there are only 40 sites





















####~~~~ Sufficentcy of Gauges 2 (NOT RUN) ~~~~####

# below is the code used when I was looking at the sufficency of gauges 2 sites
# regardless of the year of TP data. But filtring on 2001 above, I dont really care how many years of 
# data they have (right now at least...)

# # import gages 2 and filter to NYS TP CQ sites:
# 
# l.G2 <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm_EDITED_again.xlsx")
# df.G2<-reduce(l.G2, full_join, by = "STAID")
# df.G2<-select(df.G2, STAID | where(is.numeric))
# df.G2<-df.G2%>%filter(STAID %in% TP_sites$site_no)
# 
# # add the n_sample_rank: to do this:
# # rerun the temp dataframe (since I keep writing over it in different sections:
# 
# temp<-df.NWIS.TP%>%group_by(site_no)%>%
#   summarise(n=n())%>%
#   arrange(desc(n))%>%
#   mutate(n_sample_rank=rank(-n, ties.method='first'))
# 
# # left join to add n_sample rank column:
# 
# df.G2<-left_join(df.G2, temp, by = c('STAID'='site_no'))
# 
# # I want look at the number of years of data for each of these:
# 
# G2_n_years<-filter(df.NWIS.TP_CQ, site_no %in% df.G2$STAID)%>%
#   group_by(site_no)%>%
#   summarize(min_date=min(year(sample_dt)),
#             max_date=max(year(sample_dt)),
#             n_years=max(year(sample_dt))-min(year(sample_dt)))%>%
#   arrange(desc(n_years))
# 
# # I want to look at the time series of the sites with little years of data: to do this:
# 
# # filter to sits with lessthan 10 years of data:
# 
# x<-filter(G2_n_years, n_years<10)
# 
# # plot:
# 
# p<-filter(df.NWIS.TP_CQ, site_no %in% x$site_no)%>%
#   # filter(., n_sample_rank==196)%>%
#   ggplot(., aes(x = as.Date(sample_dt), y = result_va))+
#   geom_point()+
#   facet_wrap('n_sample_rank', scales ='free')
# 
# # from this plot the following sites I say have good data even though they have less than 10 years of samples:
# 
# below_10_keep<-c(46,70,72,118,119,120,122,125,126,134)
# 
# # now come up with a final list of sites: to do this:
# 
# # merge the n_samplerank and n_years dfs and filter based on n_years and the numbers of n_sample_ranks in below_10_keep:
# 
# G2_keep<-left_join(G2_n_years, temp, by = 'site_no')%>%
#   filter(., n_years>=10 | n_sample_rank %in% below_10_keep)
# 
# # make a plot:
# 
# p1<-filter(df.NWIS.TP_CQ, site_no %in% G2_keep$site_no)%>%
#   # filter(., n_sample_rank==196)%>%
#   ggplot(., aes(x = sample_dt, y = result_va))+
#   facet_wrap('n_sample_rank', scales ='free')+
#   geom_point()
# 
# # now filter the gauges2 to these sites:
# 
# df.G2<-filter(df.G2, STAID %in% G2_keep$site_no)























###############################################
 ####~~~~ Note: df.G2 has the 40 sites ~~~####
###############################################






















#####################################################
####~~~~~~~~~~~~~~~ Correlations ~~~~~~~~~~~~~~~~####
#####################################################

























#####################################################
#### Correlaitons ####
#### Using gauges2 sites ####
#####################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#### !! Editing Gauges 2 attributes !! ####

# save gauges2 df so dont need to rerun l.G2:

# df.G2.placeholder<-df.G2

# df.G2<-df.G2.placeholder

# second pass: The follow predictors have weird univariate plots,
# so removing them:

df.G2<-df.G2%>%select(-c(CONTACT, MAINS100_43, NUTR_BAS_PCT, POWER_SUM_MW, RIP100_90, RIP800_90, WOODYWETNLCD06, RAW_DIS_NEAREST_DAM))

# third pass: just because I am currious I want to remove the attributes 
# for main stem and riprian:

# df.G2<-df.G2%>%select(-starts_with('MAINS'))%>%select(-starts_with('RIP'))

# third pass not ideal

#### !! Editing Gauges 2 attributes !! ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


# combined the predictors df with the cq data:

temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.G2$STAID)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%
  left_join(., df.G2, by = c('Name'='STAID'))%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))

# create a dataframe with OLS and sens slope intercept and slopes:

df.OLS<-temp%>%
  group_by(Name)%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.I = OLS.co[1], 
            OLS.S = OLS.co[2])
  }) %>%
  ungroup

df.Sens<-temp%>%
  group_by(Name)%>%
  do({ Sens.co<-zyp.sen(log_C~log_Q,.)
  summarize(., Sen.I = Sens.co$coefficients[[1]],
            Sen.S= Sens.co$coefficients[[2]])
  }) %>%
  ungroup

# merge OLS and Sens:

df.OLS_Sens<-left_join(df.OLS, df.Sens, by = 'Name')

# merge back the watershed characteristic data to this dataframe:

df.OLS_Sens<-left_join(df.OLS_Sens, df.G2, by = c('Name'='STAID'))

# now run correlations between intercepts and slopes and watershed characteristics. to do this: (I orginally did this workflow using n_months (C:\PhD\Research\Mohawk\Code\Mohawk_Regression-analyizing_predictor_df.R)

# set up variablefor number of sites:

n_sites<-dim(df.OLS_Sens)[1] 

# use the corrr package to correlate() and focus() on your variable of choice

df.cor <- df.OLS_Sens %>% 
  # correlate() %>%
  correlate(method = 'spearman') %>%
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S))%>%
  pivot_longer(cols= c(2:5), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero

# then plot results: todo this:

# create a list of each CQ parameter (4: OLS and Sens slope and intercept)and format it for ggplotting:

l.cor<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(Spearman_Correlation)%>%  
           slice_head(n = 7)%>%
           bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           filter(!between(Spearman_Correlation, -0.25,.25))) # Order by correlation strength

# make plot list using lapply:

plist<-lapply(l.cor, \(i) i%>%ggplot(aes(x = term, y = Spearman_Correlation, color = sig_0.05)) +
                geom_bar(stat = "identity") +
                scale_color_manual(values = c("not" = "red", "sig" = "blue"),na.value = NA)+
                facet_wrap('CQ_Parameter')+
                ylab(paste('Spearman Correlation')) +
                xlab("Watershed Attribute")+
                theme(axis.text.x=element_text(angle=40,hjust=1))+
                theme(legend.position="bottom"))

# arrange plots:

p1<-ggpubr::ggarrange(plotlist = plist, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

# add title:

p1<-annotate_figure(p1, top = text_grob("Gauges 2", 
                                      color = "red", face = "bold", size = 14))
# plot:

p1

# save(p1, file = 'Processed_Data/p1.Rdata')

# univariate plots of the top correlates/small land use types:

# Lets do plots of OLS slope: to do this:

# determine the top correlates (the numberof the list determined which CQ parameter)
# here 1 = OLS intercept:

OLS<-l.cor[[1]]%>%arrange(desc(Spearman_Correlation))

# make univariate plots (facets) of different land uses and OLS slopes:
# note the predictor column isturned into an ordered factor based on thespearman correlation values
# to makethe facet plots in order:

df.OLS%>%left_join(.,df.G2%>%select(c(STAID, OLS$term)), by = c('Name'='STAID'))%>%
  pivot_longer(cols = c(4:last_col()), names_to = 'Type', values_to = 'Value')%>%
  drop_na(Value)%>%
  # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
  left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
  mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
  ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
  geom_smooth(method = 'lm')+
  geom_point()+
  facet_wrap('Type', scales = 'free')+
  ggtitle(paste('Gauges 2:', OLS$CQ_Parameter[1]))

# create a function to use lapply to plot all 4 univariate sets:

number<-3

make_plot<-function(number){
  
  OLS<-l.cor[[number]]%>%arrange(desc(Spearman_Correlation))
  
  x<-df.OLS_Sens[,c(1:5)]%>%left_join(.,df.G2%>%select(c(STAID, OLS$term)), by = c('Name'='STAID'))%>%
    pivot_longer(cols = c(6:last_col()), names_to = 'Type', values_to = 'Value')%>%
    drop_na(Value)%>%
    # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
    left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
    mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))
  
  ggplot(x, aes(x = Value, y = !!sym(OLS$CQ_Parameter[number])))+
    geom_smooth(method = 'lm')+
    geom_point()+
    facet_wrap('Type', scales = 'free')+
    ggtitle(paste('Gauges 2:', OLS$CQ_Parameter[number]))
  
}


lapply(c(1:4), \(i) make_plot(i))

#





















#####################################################
#### Correlaitons ####
#### Using Data Layers ####
#####################################################

# build a dataframe of datalayers predictors(NLCD 2001, DEM, CAFO):
# note thename df.datalayers was used above to subset the df.sf.NWIS
# to the sites after 2001. But gauges 2 cameintoplay... just usethis:

df.datalayers<-left_join(df.NWIS.NLCD.2001,df.NWIS.DEM, by = 'Name')%>%
  # left_join(.,df.NWIS.Climate, by = 'Name')%>%
  # left_join(.,df.sf.NWIS.keep, by = 'Name')%>%
  ungroup()%>%
  filter(Name %in% df.G2$STAID)%>%
  left_join(., df.sf.NWIS.keep%>%select(Name, CAFO_count), by = 'Name')%>%
  select(-geometry)

# set NA to zero:

df.datalayers[is.na(df.datalayers)]<-0

# make sure land use adds up to 100%:

sort(rowSums(df.datalayers[,-c(1, 9:12)], na.rm = T))

# looks good

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#### !! Editing data layers sites !! ####

# I want to try removing sites with over 20% developed:

# df.datalayers<-filter(df.datalayers, Developed <.2)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


# now run through correlation workflow:

temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.datalayers$Name)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%
  left_join(., df.datalayers, by = 'Name')%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))
df.OLS<-temp%>%
  group_by(Name)%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.I = OLS.co[1], 
            OLS.S = OLS.co[2])
  }) %>%
  ungroup
df.Sens<-temp%>%
  group_by(Name)%>%
  do({ Sens.co<-zyp.sen(log_C~log_Q,.)
  summarize(., Sen.I = Sens.co$coefficients[[1]],
            Sen.S= Sens.co$coefficients[[2]])
  }) %>%
  ungroup
df.OLS_Sens<-left_join(df.OLS, df.Sens, by = 'Name')
df.OLS_Sens<-left_join(df.OLS_Sens, df.datalayers, by = 'Name')
n_sites<-dim(df.OLS_Sens)[1] 
df.cor <- df.OLS_Sens %>% 
  # correlate() %>%
  correlate(method = 'spearman') %>%
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S))%>%
  pivot_longer(cols= c(2:5), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero
l.cor<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(Spearman_Correlation)%>%  
           slice_head(n = 7)%>%
           bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           distinct(term, .keep_all = T)%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           filter(!between(Spearman_Correlation, -0.25,.25))) # Order by correlation strength
plist<-lapply(l.cor, \(i) i%>%ggplot(aes(x = term, y = Spearman_Correlation, color = sig_0.05)) +
                geom_bar(stat = "identity") +
                scale_color_manual(values = c("not" = "red", "sig" = "blue"),na.value = NA)+
                facet_wrap('CQ_Parameter')+
                ylab(paste('Spearman Correlation')) +
                xlab("Watershed Attribute")+
                theme(axis.text.x=element_text(angle=40,hjust=1))+
                theme(legend.position="bottom"))
p1<-ggpubr::ggarrange(plotlist = plist, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
p1<-annotate_figure(p1, top = text_grob("Data Layers", 
                                        color = "red", face = "bold", size = 14))

p1

# multivariate plots:

OLS<-l.cor[[1]]%>%arrange(desc(Spearman_Correlation))
df.OLS_Sens%>%
  pivot_longer(cols = c(6:last_col()), names_to = 'Type', values_to = 'Value')%>%
  drop_na(Value)%>%
  select(c(Name, OLS$CQ_Parameter[1],Type, Value))%>%
  filter(Type %in% OLS$term)%>%
  # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
  left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
  mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
  ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
  geom_smooth(method = 'lm')+
  geom_point()+
  facet_wrap('Type', scales = 'free')+
  ggtitle(paste('Data Layers:', OLS$CQ_Parameter[1]))

#
number<-4

make_plot<-function(number){
  
  OLS<-l.cor[[number]]%>%arrange(desc(Spearman_Correlation))
  
  x<-df.OLS_Sens%>%
    pivot_longer(cols = c(6:last_col()), names_to = 'Type', values_to = 'Value')%>%
    drop_na(Value)%>%
    select(c(Name, OLS$CQ_Parameter[1],Type, Value))%>%
    filter(Type %in% OLS$term)%>%
    # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
    left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
    mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))
  
  ggplot(x, aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
    geom_smooth(method = 'lm')+
    geom_point()+
    facet_wrap('Type', scales = 'free')+
    ggtitle(paste('Data Layers:', OLS$CQ_Parameter[1]))
}


lapply(c(1:4), \(i) make_plot(i))

# interestingly from these univariate plots,
# it looks like if I removed sites with over 20% developed
# that would be an insane predictor of intercepts:

# I tried this and it didnt work. (It also reduced the number of sites)


























#####################################################
#### Correlaitons ####
#### Combining Gauges2 and Data Layers ####
#####################################################

# Combine the data layers and gauges 2 dfs:

df.G2_dl<-left_join(df.datalayers, df.G2, by = c('Name'='STAID'))

# looks good

# now run through correlation workflow:

temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.G2_dl$Name)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%
  left_join(., df.G2_dl, by = 'Name')%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))
df.OLS<-temp%>%
  group_by(Name)%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.I = OLS.co[1], 
            OLS.S = OLS.co[2])
  }) %>%
  ungroup
df.Sens<-temp%>%
  group_by(Name)%>%
  do({ Sens.co<-zyp.sen(log_C~log_Q,.)
  summarize(., Sen.I = Sens.co$coefficients[[1]],
            Sen.S= Sens.co$coefficients[[2]])
  }) %>%
  ungroup
df.OLS_Sens<-left_join(df.OLS, df.Sens, by = 'Name')
df.OLS_Sens<-left_join(df.OLS_Sens, df.G2_dl, by = 'Name')
n_sites<-dim(df.OLS_Sens)[1] 
df.cor <- df.OLS_Sens %>% 
  # correlate() %>%
  correlate(method = 'spearman') %>%
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S))%>%
  pivot_longer(cols= c(2:5), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero
l.cor<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(Spearman_Correlation)%>%  
           slice_head(n = 7)%>%
           bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           distinct(term, .keep_all = T)%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           filter(!between(Spearman_Correlation, -0.25,.25))) # Order by correlation strength


# make plot list using lapply:

plist<-lapply(l.cor, \(i) i%>%ggplot(aes(x = term, y = Spearman_Correlation, color = sig_0.05)) +
                geom_bar(stat = "identity") +
                scale_color_manual(values = c("not" = "red", "sig" = "blue"),na.value = NA)+
                facet_wrap('CQ_Parameter')+
                ylab(paste('Spearman Correlation')) +
                xlab("Watershed Attribute")+
                theme(axis.text.x=element_text(angle=40,hjust=1))+
                theme(legend.position="bottom"))

# arrange plots:

p1<-ggpubr::ggarrange(plotlist = plist, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
p1<-annotate_figure(p1, top = text_grob("Gauges 2 and Data Layers", 
                                        color = "red", face = "bold", size = 14))

p1

# save(p1, file = 'Processed_Data/p1.Rdata')

# univariate plots:

OLS<-l.cor[[1]]%>%arrange(desc(Spearman_Correlation))
df.OLS_Sens%>%
  pivot_longer(cols = c(6:last_col()), names_to = 'Type', values_to = 'Value')%>%
  drop_na(Value)%>%
  select(c(Name, OLS$CQ_Parameter[1],Type, Value))%>%
  filter(Type %in% OLS$term)%>%
  # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
  left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
  mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
  ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
  geom_smooth(method = 'lm')+
  geom_point()+
  facet_wrap('Type', scales = 'free')+
  ggtitle(paste('Gauges 2 and Data Layers:', OLS$CQ_Parameter[1]))

#

number<-3

make_plot<-function(number){
  
  OLS<-l.cor[[number]]%>%arrange(desc(Spearman_Correlation))
  
  x<-df.OLS_Sens%>%
    pivot_longer(cols = c(6:last_col()), names_to = 'Type', values_to = 'Value')%>%
    drop_na(Value)%>%
    select(c(Name, OLS$CQ_Parameter[1],Type, Value))%>%
    filter(Type %in% OLS$term)%>%
    # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
    left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
    mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))
  
  ggplot(x, aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
    geom_smooth(method = 'lm')+
    geom_point()+
    facet_wrap('Type', scales = 'free')+
    ggtitle(paste('Gauges 2 (edited) + Data Layers:', OLS$CQ_Parameter[1]))
}

lapply(c(1:4), \(i) make_plot(i))

#






























#######################
 ####~~~~ MLR ~~~~####
#######################

save(df.cor, df.OLS_Sens, file = 'Processed_Data/troubleshoot_AIC.Rdata')

# create the l..cor list as above but now only filter out the 
# non-significant attributes:

load('Processed_Data/troubleshoot_AIC.Rdata')

l.cor.MLR<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(Spearman_Correlation)%>%  
           # slice_head(n = 7)%>%
           # bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           # distinct(term, .keep_all = T)%>%
           # filter(!between(Spearman_Correlation, -0.25,.25))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           # mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           # filter(sig_0.05 == 'sig')%>%
           as.data.frame())


x<-l.cor.MLR[[1]]

# create a list of dataframes from df.OLS.Sens and subseting the attributes using the 
# names in each l.cor[[i]]$term for simple and full models:

l.cor.MLR.full<-lapply(1:4, \(i) df.OLS_Sens%>%
                    select(i+1, l.cor.MLR[[i]]$term)%>%
                    as.data.frame()%>%
                      rename(term = 1))

x<-l.cor.MLR.full[[1]]

# loop through data frames and make lm and stepAIC objects (wont work on lapply for some reason):

# set up list for aic objects to append into:

l.aic<-list()

# loop:

i<-1

for (i in 1:4){
  
  # find the best predictor name for the simple model:
  
  n<-l.cor.MLR[[i]]$term[which.max(abs(l.cor.MLR[[i]]$Spearman_Correlation))]
  
  # createformula with this predictor:
  
  n<-paste('term ~', n)
  
  # create simple model:
  
  m.simple<-lm(n,l.cor.MLR.full[[i]])
  
  # create full model:
  
  m.full<-lm(term~., l.cor.MLR.full[[i]])
  
  # create AIC object:
  
  aic<-stepAIC(m.simple, scope = list(upper=m.full, lower =~1), direction = 'both', trace = FALSE)

  # save aic object to list:
  
  l.aic[[i]]<-aic
  
}

# look at summary of aic objects:

lapply(l.aic, summary)

# set names of list elements:

names(l.aic)<-names(df.OLS_Sens)[2:5]

names(l.aic[[1]]$model)[-1]
  
# extract best predictors for each parameter:

l<-lapply(l.aic, \(i) names(i$model)[-1])

# look at univariate plots of these. to do this:

# make a list of ggplot objects for each parameter:

l.MLR.plots<-lapply(1:4, \(i) df.OLS_Sens%>%
                       select(i+1, l[[i]])%>%
                       as.data.frame()%>%
                      pivot_longer(cols = 2:last_col(), names_to = 'Attribute', values_to = 'value')%>%
                      mutate(Attribute = factor(Attribute, levels = l[[i]]))%>%
                      ggplot(., aes(y = !!sym(names(df.OLS_Sens)[i+1]), x = value))+
                      facet_wrap('Attribute', scales = 'free')+
                      geom_point()+
                      geom_smooth(method = 'lm')+
                      ggtitle((names(df.OLS_Sens)[i+1]))
                       
                    
                    )

l.MLR.plots

#



















#### Categorizing land use ####

# USGS criteria:
# Agricultural sites have >50% agricultural land and ≤5% urban land;
# urban sites have >25% urban and ≤25% agricultural land; 
# undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land; 
# all other combinations of urban, agricultural, and undeveloped lands are classified as mixed

# combine Ag and Pasture into a single landuse for Ag:

df.datalayers<-mutate(df.datalayers, Ag2 = Ag+Pasture)

# set NA to zero

df.datalayers[is.na(df.datalayers)]<-0

# create the land use class column based on USGS critiera:
# adjusted thresholds:

df.datalayers<-df.datalayers%>%
  mutate(USGS.LU.Adjusted = 'Mixed')%>%
  mutate(USGS.LU.Adjusted = case_when(.default = 'Mixed',
                                      Ag2 > .30 & Developed <= .1 ~ 'Agriculture',
                                      Developed > .1 & Ag2 <= .3 ~ 'Urban',
                                      Developed <= .1 & Ag2 <= .1 ~ 'Undeveloped'))

# merge the OLS and Sens slopes and intercepts with this df:

df.datalayers<-left_join(df.datalayers, df.OLS_Sens[,1:5], by = 'Name')

# plot

p<-df.datalayers%>%
  pivot_longer(cols = 15:18, names_to = 'CQ_parameter', values_to = 'Value')%>%
  ggplot(., aes(x=USGS.LU.Adjusted, y=Value, color =USGS.LU.Adjusted ))+
    geom_boxplot(varwidth = TRUE, alpha=0.2)+
    # scale_x_discrete(labels=my_xlab)+
    facet_wrap('CQ_parameter', scales = 'free')+
    stat_compare_means(method = "anova", label.y = max(df.datalayers$Value))+      # Add global p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = "0.5") +
    ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')

p

#






















#### Categorizing land use (Not run) - experimenting with different data layers ####
# 
# # USGS criteria:
# # Agricultural sites have >50% agricultural land and ≤5% urban land;
# # urban sites have >25% urban and ≤25% agricultural land; 
# # undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land; 
# # all other combinations of urban, agricultural, and undeveloped lands are classified as mixed
# 
# # I have four different datasets to pull land use from:
# # NLCD 2001 and 2019
# # CDL 2008 and 2020
# 
# # The plan is to make two facet plots, one for orginal USGS thresholds and one for adjusted USGS thresholds:
# # facets will be the different CQ parameters with x axis being the different land use catgeory andyaxis beung the parameter value
# # for each land use, color will be used to determine which dataset the values came from
# 
# 
# # In a first pass using these thresholds the number of ag and urban sites wasvery low
# # I will play with these numbers to see what happens
# 
# # merge land use from CDL and NLCD to df.sf.NWIS.keep.2 and add a yearcolumn:
# 
# t1<-left_join(df.sf.NWIS.keep.2, df.NWIS.CDL%>%mutate(LU_source = 'CDL_2020'), by = 'Name')
# t2<-left_join(df.sf.NWIS.keep.2, df.NWIS.CDL.2008%>%mutate(LU_source = 'CDL_2008'), by = 'Name')
# t3<-left_join(df.sf.NWIS.keep.2, l.NWIS.NLCD$`2019`%>%mutate(LU_source = 'NLCD_2019'), by = 'Name')
# t4<-left_join(df.sf.NWIS.keep.2, l.NWIS.NLCD$`2001`%>%mutate(LU_source = 'NLCD_2001'), by = 'Name')
# 
# # merge the CDL and NLCD dfs:
# 
# df.LU<-bind_rows(t1,t2,t3,t4)
# 
# # combine Ag and Pasture into a single landuse for Ag:
# 
# df.LU<-mutate(df.LU, Ag = Ag+Pasture)
# 
# # set NA to zero
# 
# df.LU[is.na(df.LU)]<-0
# 
# # create the land use class column based on USGS critiera:
# 
# # orginal thresholds:
# 
# df.LU<-df.LU%>%
#   mutate(USGS.LU = 'Mixed')%>%
#   mutate(USGS.LU = case_when(.default = 'Mixed',
#                              Ag > .50 & Developed <= .05 ~ 'Agriculture',
#                              Developed > .25 & Ag <= .25 ~ 'Urban',
#                              Developed <= .05 & Ag <= .25 ~ 'Undeveloped')
#   )
# 
# # adjusted thresholds:
# 
# df.LU<-df.LU%>%
#   mutate(USGS.LU.Adjusted = 'Mixed')%>%
#   mutate(USGS.LU.Adjusted = case_when(.default = 'Mixed',
#                              Ag > .30 & Developed <= .1 ~ 'Agriculture',
#                              Developed > .1 & Ag <= .3 ~ 'Urban',
#                              Developed <= .1 & Ag <= .1 ~ 'Undeveloped'))
# 
# # merge the OLS and Sens slopes and intercepts with this df:
# 
# df.LU<-left_join(df.LU, df.OLS_Sens[,1:5], by = 'Name')
# 
# # pivot longer for geom_box + facet:
# 
# data<-df.LU%>%
#   pivot_longer(cols = 16:19, names_to = 'CQ_parameter', values_to = 'Value')%>%
#   mutate(USGS.LU=factor(USGS.LU))
# 
# # make ggplot:
# 
# # for orginal thresholds:
# 
# my_xlab <- paste(levels(factor(df.LU$USGS.LU)),"\n(N=",table(factor(df.LU$USGS.LU)),")",sep="")
# 
# ggplot(data, aes(x=LU_source, y=Value, color =USGS.LU ))+
#   geom_boxplot(varwidth = TRUE, alpha=0.2)+
#   # scale_x_discrete(labels=my_xlab)+
#   facet_wrap('CQ_parameter', scales = 'free')+
#   stat_compare_means(method = "anova", label.y = max(data$Value))+      # Add global p-value
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "0.5") +
#   ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')
# 
# # for adjusted thresholds:
# 
# my_xlab <- paste(levels(factor(df.LU$USGS.LU.Adjusted)),"\n(N=",table(factor(df.LU$USGS.LU.Adjusted)),")",sep="")
# 
# ggplot(data, aes(x=LU_source, y=Value, color =USGS.LU.Adjusted ))+
#   geom_boxplot(varwidth = TRUE, alpha=0.2)+
#   # scale_x_discrete(labels=my_xlab)+
#   facet_wrap('CQ_parameter', scales = 'free')+
#   stat_compare_means(method = "anova", label.y = max(data$Value))+      # Add global p-value
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "0.5") +
#   ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')
# 

























#### Grouping CQ curves (stationary, mobilization, dilutionary, complex) ####

# create a list of dataframes for each sites CQ observations:

temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.datalayers$Name)%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))

l.temp<-temp%>%
  split(., .$Name)

# create lm models for each site:

l.lm.CQ_slopes<-lapply(l.temp, \(i) lm(log_C~log_Q, data=i))

# save the model coef ad pvals:

coef<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,1] ))), 'site_no')

# save the pvalues 

pvals<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,4] ))), 'site_no')%>%
  rename(I.pval = 2, S.pval = 3)

# merge the two dfs:

m<-left_join(coef,pvals,by='site_no')

# add column for CQ type:

m<-mutate(m, Type = ifelse(S.pval>0.05, 'Stationary', ifelse(log_Q>0, 'Mobilization', 'Dilution')))

# merge labels with plotting df:

temp<-left_join(temp,m%>%select(site_no, Type),by=c('Name'='site_no'))

# create a df for CQplot of all sites with BP analysis:

df_Seg.2<-filter(df_Seg, site %in% temp$Name)%>%
 left_join(.,m%>%select(site_no, Type),by=c('site'='site_no')) 

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'purple', size = 2)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# looking at this plot I want to add a fourth CQ type for complex, if the slopes of the BP analysis look widely different. 
# I will start wit hcalcuating the angle between pre-post BP slope and see which sites get signled out:

# add a new column with the angle between the two lines:

df_Seg.2<-df_Seg.2%>%
  mutate(slope_angle=factor(round(atan(abs((Slope2-Slope1)/(1+(Slope2*Slope1)))),1)))

# create color pallete for the slope angle:

hc<-heat.colors(length(unique(df_Seg.2$slope_angle)), rev = T)

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  new_scale_color() +
  geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
  geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
  scale_color_manual(name = "Slope Angle", values = hc)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

p

# based on this plot, I would chose the following sites as complex:

complex_sites<-unique(df_Seg.2$site)[c(3,16,21,22,29,37)] 

# old complex site (when n=53): c(3,18,24,25,26,35,42,43)

df_Seg.2<-mutate(df_Seg.2, Type = ifelse(site %in% complex_sites, 'Complex', Type))

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# I want to add land use as color add number of CAFOs as text. to do this:

# merge the plotting df with the land use for NLCD 2001 and CAFO count and adjusted thresholds:

df_Seg.2<-left_join(df_Seg.2, df.datalayers%>%select(Name, USGS.LU.Adjusted, CAFO_count), by = c('site'='Name'))

# if CAFO count is zero set to NA:

df_Seg.2$CAFO_count[df_Seg.2$CAFO_count==0]<-NA

# plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

# save(p1, file = 'Processed_Data/p1.Rdata') 
  
# add CQ type to the mapping df with polygon trasperncy based on number of samples:

map.CQ_Type<-df.sf.NWIS.keep%>%
  filter(Name %in% df.datalayers$Name)%>%
  left_join(.,distinct(df_Seg.2, site, .keep_all = T)%>%select(.,c(site, Type, n_sample_rank)),  by = c('Name'='site'))%>%
  select(Name, Type, n_sample_rank)%>%
  arrange(n_sample_rank)%>%
  mutate(NEW = (1/row_number())*40)
  # left_join(.,df.NWIS.TP%>%group_by(site_no)%>%summarise(n=n()), by = c('site'='site_no'))%>%
  # mutate(NEW = case_when(n < 50 ~ .05,
  #                        n >=50 & n < 100 ~ .25,
  #                        n >=100 & n < 500 ~ .5,
  #                        n >=500 ~ .9))
# map:

# mapview(map.CQ_Type, zcol = 'Type', alpha.regions = 'NEW')

######################################################
 #### ~~~~ XY plot of Slopes and Intercepts ~~~~ ####
####################################################

# plot of slopes and intercepts for each CQ type:

m<-m%>%mutate(Type = ifelse(site_no %in% complex_sites, 'Complex', Type))%>%
  left_join(., df.datalayers%>%select(Name, USGS.LU.Adjusted), by = c('site_no'='Name'))

ggplot(m, aes(x=`(Intercept)`, y=log_Q, color = Type))+
  geom_point(size = 2)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  facet_wrap('USGS.LU.Adjusted', scales = 'fixed')
  
# I want to do this with more sites























#### Land Use changes across watersheds (not run) ####

# # Look at changes in NLCD and CDL calculated Ag %land between 2008-2020 for each 53 watershed
# # to do this:
# 
# # pivot longerto getthe land use catgories in a single column:
# 
# df.LU.2<-pivot_longer(df.LU, cols = 6:12, names_to = 'LandUse', values_to = 'values')
# 
# # pivot wider to get the land use source (data layers) in their own columns:
# 
# df.LU.2<-pivot_wider(df.LU.2, names_from = 'LU_source', values_from = 'values')
# 
# # pivot longer to put CDL and NLCD in own columns:
# 
# df.LU.CDL<-pivot_longer(df.LU.2, cols = starts_with('CDL'), names_to = 'CDL', values_to = 'CDL_values')
# df.LU.NLCD<-pivot_longer(df.LU.2, cols = starts_with('NLCD'), names_to = 'NLCD', values_to = 'NLCD_values')
# 
# # group by land use and calculate the difference between the years:
# 
# df.LU.2<-df.LU.2%>%mutate(CDL_diff = abs(round(CDL_2020-CDL_2008,2)),NLCD_diff = abs(round(NLCD_2019-NLCD_2001,2)))%>%
#   arrange(desc(CDL_diff))
# 
























#### Seasonal ANanlysis ####

# is the CQ relaitonship different between seasons at these sites?

# I need to filter down to sites with over 100 samples for this, 20 samples isnt going tocut it:

df_Seg.3<-filter(df_Seg.2, n>100)

# add season column:

df_Seg.3$Season<-getSeason(df_Seg.3$Date)

# plot:

p1<-ggplot(df_Seg.3, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Season), size = 4)+
  # scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.title=element_text(size=14)
  )+
  geom_rect(data = df_Seg.3%>%distinct(df_Seg.3$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = Type), alpha = .15)+
  scale_fill_manual(name = "CQ Type", values = c("red", "blue","purple", "green"))

p1

# save(p1, file = 'Processed_Data/p1.Rdata') 

# 


























############################################
 #### ~~~~ Q Exceedence Proability ~~~~ ####
############################################




























####~~~~ Calcualte EP ~~~~####

# subset the daily flows to the 40 sites:

temp<-filter(df.NWIS.Q.for_TN, site_no %in% df.datalayers$Name)

# split into lists by site:

temp<-split(temp, f = temp$site_no) 

# sort each dataframe in the list by the flow, add a column for m and n, convert back to a dataframe using bind_rows, and select only the needed columns for the left join (next step):

temp<-bind_rows(lapply(temp, \(i) i[order(i$X_00060_00003,decreasing = T),]%>%mutate(m = 1:n(), n = n())))%>%select(site_no, Date, m, n)

# append the value of M and n for each C-Q observation in the C-Q dataframe:
# need to rename the n column in the left dataframe as well:

df_Seg.2<-left_join(df_Seg.2%>%rename(n_samples = n), temp, by = c("site" = "site_no", "Date" = "Date"))

# calcualte the exceednce proability for Q for each C observation:

df_Seg.2$EP<-round(df_Seg.2$m/(df_Seg.2$n+1), 4)

# plot:

p1<-ggplot(df_Seg.2, aes(x = EP, y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "Log-Log\nCQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  # geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = .5, y =-2, label = CAFO_count), inherit.aes = FALSE, size=15)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  scale_x_reverse()+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p1

# save(p1, file = 'Processed_Data/p1.Rdata')


























####~~~~ Estimate Break point analysis using EP-Q ~~~~####

#  rerun the BP analysis using EP-Q:

l_Seg.EP<-list()
davies.test.matrix.EP<-NULL
EP<-df.NWIS.TP_CQ%>%
  filter(site_no %in% df.datalayers$Name)%>%
  left_join(., df_Seg.2%>%select(site, Date, EP), by = c('site_no'='site', 'sample_dt'='Date'))%>%
  drop_na(result_va, Q_yield)%>%
  rename(site=site_no, Date = sample_dt)
temp.loop<-sort(unique(EP$site))

i<-1

for (i in seq_along(temp.loop)){
  
  tryCatch({
    
    # print the site name for loop debugging:
    
    print(i)
    print(temp.loop[i])
    
    # create a dataframe that will work with segmented. To do this: 
    # filter for the site the loop is in
    # add log transformed C and Q columns, as well as duplicated columns for renamed C and Q
    # filter for real log C and Q values so breakpoint analysis works smoothly:
    
    df<-EP%>%
      filter(site == temp.loop[i])%>%
      mutate(C = result_va)%>%
      mutate(log_C = log(C))%>%
      filter(is.finite(log_C))%>%
      filter(is.finite(EP))
    
    # build a single slope lm for log C and Q. Tis model is also used inthe breakpoint analysis inthenext step:
    
    m<-lm(log_C~EP, df)
    
    # perform breakpoint regression:
    
    m_seg<-segmented(obj = m, npsi = 1)
    
    # perform davies test for constant linear predictor:
    # the results are saved as a string with the the site name and true/false:
    
    x<-paste(temp.loop[i], '-', davies.test(m)$p.val<0.05)
    
    # add the results of davies test to the matrix made prior to this for loop:
    
    davies.test.matrix.EP<-c(davies.test.matrix.EP,x)
    
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
    
    
    # get the model fitted data, put in a dataframe, and reformat this dataframe to export out of the loop
    
    if(length(class(m_seg))==2){
      fit <- data.frame(site = df$site, Date = df$Date, Seg_C_EP = fitted(m_seg), I1_EP = inter$Est.[1], I2_EP = inter$Est.[2], Slope1_EP = s[1,1], Slope2_EP = s[2,1], BP_EP = bp)
      result_df<-left_join(df, fit, by = c('site', 'Date'))
    } else{
      fit <- data.frame(site = df$site, Date = df$Date, Seg_C_EP = fitted(m_seg), I1_EP = inter$Est.[1], I2_EP = NA, Slope1_EP = NA, Slope2_EP = NA, BP_EP = NA)
      result_df<-left_join(df, fit, by = c('site', 'Date'))
    }
    
    l_Seg.EP[[i]]<-result_df
    
  }, 
  
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

davies.test.matrix.EP
df.davies.EP<-as.data.frame(davies.test.matrix.EP)%>%
  separate_wider_delim(1, "-", names = c("site", "BP_EP_yes"))%>%
  mutate(across(c(1,2), trimws))
df_Seg.EP<-bind_rows(l_Seg.EP)%>%
  mutate(Seg_EP=EP)
df_Seg.EP<-df_Seg.EP%>%
  left_join(., df.davies.EP, by = 'site')%>%
  mutate(across(c(Seg_C_EP, Seg_EP), ~replace(., BP_EP_yes == 'FALSE', NA)))

# merge the land use, CQ type, data to this df:

df_Seg.EP<-left_join(df_Seg.EP, m%>%select(site_no, Type, USGS.LU.Adjusted), by = c('site'='site_no'))

# plot:

p<-ggplot(df_Seg.EP, aes(x = EP, y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "Log-log\nCQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  new_scale_color()+
  geom_line(aes(x = Seg_EP, y = Seg_C_EP), color = 'yellow', size = 1.5)+
  scale_color_manual(name = "Log(C)~EP-Q\n Breakpoint Analysis", values = "yellow")+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free_y')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  scale_x_reverse()+
  geom_rect(data = df_Seg.EP%>%distinct(df_Seg.EP$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

save(p, file = 'Processed_Data/p1.Rdata')

# I like the idea of using the EP plots to help ID which
# sites are complex. For example, site 22 was deemed complex
# with log-log plots, but in log-EP space it looks similar to
# the other mobilization sites. 























#----------------------------------------#
  ####~~~  Reclassifying CQ type ~~~#### 
   ###~~~  using the EP-Q plots  ~~~### 
#----------------------------------------#

# rerun the df_Seg.2 workflow bt=ut changing the complex_sites vector:
temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.datalayers$Name)%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))
l.temp<-temp%>%
  split(., .$Name)
l.lm.CQ_slopes<-lapply(l.temp, \(i) lm(log_C~log_Q, data=i))
coef<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,1] ))), 'site_no')
pvals<-tibble::rownames_to_column(as.data.frame(t(sapply(l.lm.CQ_slopes, \(i) summary(i)$coefficients[,4] ))), 'site_no')%>%
  rename(I.pval = 2, S.pval = 3)
m<-left_join(coef,pvals,by='site_no')
m<-mutate(m, Type = ifelse(S.pval>0.05, 'Stationary', ifelse(log_Q>0, 'Mobilization', 'Dilution')))
temp<-left_join(temp,m%>%select(site_no, Type),by=c('Name'='site_no'))
df_Seg.2<-filter(df_Seg, site %in% temp$Name)%>%
  left_join(.,m%>%select(site_no, Type),by=c('site'='site_no')) 
df_Seg.2<-df_Seg.2%>%
  mutate(slope_angle=factor(round(atan(abs((Slope2-Slope1)/(1+(Slope2*Slope1)))),1)))
hc<-heat.colors(length(unique(df_Seg.2$slope_angle)), rev = T)
complex_sites<-unique(df_Seg.2$site)[c(3,16,21,37)] 
df_Seg.2<-mutate(df_Seg.2, Type = ifelse(site %in% complex_sites, 'Complex', Type))
df_Seg.2<-left_join(df_Seg.2, df.datalayers%>%select(Name, USGS.LU.Adjusted, CAFO_count), by = c('site'='Name'))
df_Seg.2$CAFO_count[df_Seg.2$CAFO_count==0]<-NA

# plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type), size = 1.5)+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  # geom_text(aes(x = 2, y = -2, label = CAFO_count), inherit.aes = FALSE, size = 30)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .15)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","purple", "green"))

p

#

























































# finally save image to workspace:

# save.image(file = 'Processed_Data/NWIS.Rdata')

load('Processed_Data/NWIS.Rdata')




