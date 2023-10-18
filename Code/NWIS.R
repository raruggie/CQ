# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

library(climateR)
library(ggpubr)
library(ggnewscale)
library(readxl)
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

####################### Process the NWIS database for the CQ analysis #######################

# this code adapted from NYS_site_ranking.R

# Download the metadata for sites with daily flow data. To do this:

# Use dataRetrieval::whatNWISdata to query USGS gauge sites with daily flow data in NY

df.NWIS.Q_sites<- whatNWISdata(stateCd = 'NY', parameterCd = "00060")

write.csv(df.NWIS.Q_sites, "C:/PhD/CQ/Raw_Data/df.NWIS.Q_sites.csv", row.names=FALSE)

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

df.NWIS.TP_sites<-fun.df.Pair_consit_flow('00665', df.NWIS.Q_sites, n_samples = 20)

write.csv(df.NWIS.TP_sites, "C:/PhD/CQ/Raw_Data/df.NWIS.TP_sites.csv", row.names=FALSE)

# download the raw daily flow data for these sites. To do this:
# readNWISdv can be used to download a single dataframe with all the raw flow data for a vector of gauge numbers (this takes a looong time):

df.NWIS.Q<-readNWISdv(siteNumbers = df.NWIS.TP_sites$site_no, parameterCd = '00060', startDate = "", endDate = "", statCd = "00003")

write.csv(df.NWIS.Q, "C:/PhD/CQ/Raw_Data/df.NWIS.Q.csv", row.names=FALSE)

# download the raw discrete TP sample data:

df.NWIS.TP<-readNWISqw(siteNumbers = df.NWIS.TP_sites$site_no, parameterCd = '00665')

write.csv(df.NWIS.TP, "C:/PhD/CQ/Raw_Data/df.NWIS.TP.csv", row.names=FALSE)

# now join the TP with the flow df:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP, df.NWIS.Q, by=c("site_no"="site_no", "sample_dt"="Date"))

# remove observations where there are not CQ pairs:

df.NWIS.TP_CQ<-df.NWIS.TP_CQ%>%drop_na(X_00060_00003)

# arrange by number of TP observations. To do this:

# first arrange the dataframe with the number of samples: 

temp<-df.NWIS.TP%>%group_by(site_no)%>%
  summarise(n=n())%>%
  arrange(desc(n))

# add a column with the rank:

temp$n_sample_rank<-1:nrow(temp)

# finally merge this df with df.NWIS.TP_CQ and arrange by the new column:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP_CQ,temp, by='site_no')%>%
  arrange(n_sample_rank)

# next step is to use Q yield to get better looking plot... idk if it will help but want to try also doesnt hurt to have the watershed areas as well. To do this:

# download the drainage areas from site metadata using readNWISsite

df.NWIS.TP_site_metadata<-readNWISsite(siteNumbers = unique(df.NWIS.TP_CQ$site_no))

write.csv(df.NWIS.TP_site_metadata, "C:/PhD/CQ/Raw_Data/df.NWIS.TP_site_metadata.csv", row.names=FALSE)

# then select just the site number and DA column:

df.DA<-df.NWIS.TP_site_metadata%>%
  select(site_no, drain_area_va)

# finally merge with df.NWIS.TP_CQ and create a new Q column with area normalized flows (not worrying about units right now): 
# Note: will filter for NA in C and Q for breakpoint analysis in the next step as to keep the full list of sites with CQ pairs in this dataframe.
# Note: Some sites returned NA on draiange areas in readNWISsite, but I'll delinate anyways so I want the full list:

df.NWIS.TP_CQ<-left_join(df.NWIS.TP_CQ, df.DA, by = 'site_no')%>%
  mutate(Q_yield = X_00060_00003/drain_area_va)

# I want to to fit breakpoint models to the CQ curves. To do this:

# create empty list to hold the results of the segmented function (which does the breakpoint analysis)

l_Seg<-list()

# create a matrix to hold the results of the davies test, which determines if a two slope model is warrented over a single slope model:

davies.test.matrix<-NULL

# create a new dataframe of only paired CQ observations (such that the breakpoit analysis function runs smoothly) (I didnt want to do this in loop for some reason, I cant remeber why but it wouldnt work):

df.NWIS.TP_CQ_for_BP<-df.NWIS.TP_CQ%>%
  drop_na(result_va, Q_yield)

# create a vector of ordered unique site names:

temp.loop<-sort(unique(df.NWIS.TP_CQ_for_BP$site_no))

# test i for for loop building:

# i<-unique(df.NWIS.TP_CQ_for_BP$site_no)[4]

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
    
    df<-df.NWIS.TP_CQ_for_BP%>%
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
      result_df<-fit%>%mutate(site = temp.loop[i], Q_real = df$Q, C = df$C, I1 = inter$Est.[1], I2 = inter$Est.[2], Slope1 = s[1,1], Slope2 = s[2,1], BP = bp)
    } else{
      result_df<-fit%>%mutate(site = temp.loop[i], Q_real = df$Q, C = df$C, I1 = NA, I2 = NA, Slope1 = NA, Slope2 = NA, BP = NA)
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

# save(map.NWIS.TP_sites, file = 'C:/PhD/CQ/Processed_Data/map.NWIS.TP_sites.Rdata')

####################### Run through CQ metrics/watershed attribute correlation #######################

# build a matrix of the number of TP samples as a function of watershed size. To do this:

# I already have a df with paired CQ observations and DA, just need to get the distinct sites:

TP_sites<-df.NWIS.TP_CQ%>%distinct(site_no, .keep_all = T)

# set up df to populate:

m<-data.frame(Min_num_samples  = c(20,50,75,100,200), '25' = NA, '50' = NA, '100' = NA, '150'=NA, '250'=NA, '500'=NA, '1000'=NA, 'Unlimited'=NA)

# set up variables for thresholds for number of samples and DA size:

n_sam<- c(20,50,75,100,200)-1

min_DA<- c(25,50,100,150,250,500,1000)

# loop through number of samples (rows):

# i<-1

for (i in seq_along(n_sam)){
  
  temp.i<-TP_sites%>%filter(n>n_sam[i])
  
  m$Unlimited[i]<-dim(temp.i)[1]

  # j<-2
  
  # loop through the size of the DA (columns)
  for(j in  seq_along(min_DA)){
    
    temp.j<-temp.i%>%filter(drain_area_va<=min_DA[j])
    
    m[i,j+1]<-dim(temp.j)[1]
    
  }
  
}

# going to use this matrix to pick subsets of the NWIS sites and see how it affects the correlations analysis:

# correlation analysis (adapated from Skan_CQ.R):

# create a single dataframe of the CQ observations and watershed attributes for all sites. to do this:

# read in the gauges 2 database (I am forgoing the datalayers workflow for now since I may have a great set of predictors in this database):

# read in all sheets using function, remove the last element (does notcomtain useful info), and convert to a single df:

l.G2 <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm.xlsx")

l.G2[[27]]<-NULL

df.G2<-reduce(l.G2, full_join, by = "STAID")

# filter the gauges2 to the NYS TP CQ sites:

df.G2<-df.G2%>%filter(STAID %in% unique(df.NWIS.TP_CQ$site_no))

# only 89 of the orginal 137 sites are in gauges 2

# lets see what the tradeoff matrix looks like if only these 89 are included:

TP_sites<-df.NWIS.TP_CQ%>%
  filter(site_no %in% df.G2$STAID)%>%
  distinct(site_no, .keep_all = T)

m<-data.frame(Min_num_samples  = c(20,50,75,100,200), '25' = NA, '50' = NA, '100' = NA, '150'=NA, '250'=NA, '500'=NA, '1000'=NA, 'Unlimited'=NA)

n_sam<- c(20,50,75,100,200)-1

min_DA<- c(25,50,100,150,250,500,1000)

# i<-1

for (i in seq_along(n_sam)){
  
  temp.i<-TP_sites%>%filter(n>n_sam[i])
  
  m$Unlimited[i]<-dim(temp.i)[1]
  
  # j<-2
  
  for(j in  seq_along(min_DA)){
    
    temp.j<-temp.i%>%filter(drain_area_va<=min_DA[j])
    
    m[i,j+1]<-dim(temp.j)[1]
    
  }
  
}

# I cant make the call, 
# iwantto just stick withthe plan, so now going back to datalayers analysis
# so now I need to delinate the 137 NWIS sites:

l.SS_WS.NWIS<-lapply(seq_along(df.NWIS.TP_site_metadata$site_no), \(i) Delineate(df.NWIS.TP_site_metadata$dec_long_va[i], df.NWIS.TP_site_metadata$dec_lat_va[i]))

names(l.SS_WS.NWIS)<-df.NWIS.TP_site_metadata$site_no

# save(l.SS_WS.NWIS, file = 'C:/PhD/CQ/Downloaded_Data/l.SS_WS.NWIS.Rdata')

df.sf.NWIS<-fun.l.SS_WS.to.sfdf(l.SS_WS.NWIS)

# save(df.sf.NWIS, file = 'C:/PhD/CQ/Processed_Data/df.sf.NWIS.Rdata')

load('C:/PhD/CQ/Processed_Data/df.sf.NWIS.Rdata')

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

temp<-df.sf.NWIS[,c(1,107)]%>%
  filter(abs(Delination_Error)>.02)

mapview(temp)
  
# after playing with the numbers, 2% seems like a good threshold for non-delineations:
# subset df.sf to those sites that are considered to delineate correctly:

df.sf.NWIS.keep<-df.sf.NWIS%>%filter(abs(Delination_Error)<=.02)

# note I suspect that snapping the lat long to nhd would improve the delineation success notgoing to do that now

#### Data Layers ####

# now run the datalayers workflow for these sites:

#### CDL: 

# note when I ran ths CDL code block I used the 103 sites thatwere apartofthe set, this is prior to filtering based onthe CDLdate... so I am just going to filter the resulting CDL df at the end. Just note that what is saved is of the 103 sites, but that the code doesnt reflect this, i.e. I just ran the filter anddidnt keep the code because df.NWIS.TP.keep would be down to 56 if I did it right 

# download:

rast.NWIS.CDL.2020 <- GetCDLData(aoi = df.sf.NWIS.keep, year = "2020", type = "b", tol_time = 1000)

# convert to rast:

rast.NWIS.CDL.2020<-rast(rast.NWIS.CDL.2020)

# plot

plot(rast.NWIS.CDL.2020)

# reproject to sample watershed vector data to match raster data:

vect.NWIS<-vect(df.sf.NWIS.keep) # need to first ininalize vect since sometimes reading in rdata file (terra issue)

vect.NWIS.proj<-terra::project(vect.NWIS, crs(rast.NWIS.CDL.2020))

# extract frequency tables for each sample watershed

l.NWIS.CDL <- terra::extract(rast.NWIS.CDL.2020, vect.NWIS.proj, table,ID=FALSE)[[1]]

# save(l.NWIS.CDL, file='Processed_Data/l.NWIS.CDL.Rdata')

load('Processed_Data/l.NWIS.CDL.Rdata')

# aggregate CDL: to do this:
# looking at the CDL legend (linkdata), I determined the following:

Ag<-c(1:6,10:14,21:39,41:61,66:72,74:77,204:214,216:227,229:250,254)
Pasture<-c(176)
Forest<-c(63,141:143)
Shrub<-c(64,152)
Barren<-c(65,131)
Developed<-c(82,121:124)
Water<-c(83,111)
Wetlands_all<-c(87,190,195)
Other<-c(88,112)

l <- tibble::lst(Ag,Pasture,Forest,Shrub,Barren,Developed,Water,Wetlands_all,Other)

reclass_CDL<-data.frame(lapply(l, `length<-`, max(lengths(l))))%>%
  pivot_longer(cols = everything(), values_to = 'MasterCat',names_to = 'Crop')%>%
  drop_na(MasterCat)

# convert resulting list of tables to list of dfs

l.NWIS.CDL<-lapply(l.NWIS.CDL, as.data.frame)

# left join each df in the list to the CDL legend key, as well as calcuate the pland:
  
l.NWIS.CDL<-lapply(l.NWIS.CDL, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., reclass_CDL, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes

# set names of list:

names(l.NWIS.CDL)<-df.sf.NWIS.keep$Name # note doesnt work with the workflow set up to filter on 2008 limiter

# combine list into single df:

df.NWIS.CDL<-bind_rows(l.NWIS.CDL, .id = 'Name')

# remove potential for one of the MasterCats in the CDL to be empty, which is messing with the pivot_wider below:

df.NWIS.CDL<-filter(df.NWIS.CDL, Crop != '')

# pivot wider:

# if using the CDL linkdata, use this:

# df.NWIS.CDL<- pivot_wider(df.NWIS.CDL[,-2], names_from = Crop, values_from = Freq)

# if using the reclass_CDL data, use this:

df.NWIS.CDL<-df.NWIS.CDL[,-2]%>%
  group_by(Name, Crop)%>%
  summarise(Freq=sum(Freq, na.rm = T))%>%
  pivot_wider(., names_from = Crop, values_from = Freq)

# check to see if add up to 100%:

sort(rowSums(df.NWIS.CDL[,-1], na.rm = T))

# looks good. 

# Done with CDL

#### NED

# download:

DEM.NWIS<-get_ned(df.sf.NWIS.keep, label = '2') # already SpatRaster!

writeRaster(DEM.NWIS, file='Downloaded_Data/DEM.NWIS.tif')

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

vect.NWIS<-vect(df.sf.NWIS.keep)

vect.NWIS.proj<-terra::project(vect.NWIS, crs(DEM.NWIS))

# extract the metrics over each watershed using the function above:

df.NWIS.DEM <- as.data.frame(terra::extract(DEM.NWIS, vect.NWIS.proj, f))

# set the names of the df:

names(df.NWIS.DEM)<-c('Name', 'Elev_Avg', 'Elev_Range', 'Elev_SD')

# set the names of the sites:

df.NWIS.DEM$Name<-df.sf.NWIS.keep$Name

# finally save the df:

save(df.NWIS.DEM, file= 'Processed_Data/df.NWIS.DEM.Rdata')

# done with DEM

#### Climate

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

## combine climate,CDL, DEM, and SS_WS to one df:

df.NWIS.Predictors<-left_join(df.NWIS.CDL,df.NWIS.DEM, by = 'Name')%>%left_join(.,df.NWIS.Climate, by = 'Name')%>%left_join(.,df.sf.NWIS.keep, by = 'Name')

# of these sites, Iwant to filter down those that havedata during the CDLperiod
# the CDL goes back to 2008: todothis:

# filter theraw TPdatabased on date and toget the desired sites:

temp<-df.NWIS.TP_CQ%>%filter(year(df.NWIS.TP_CQ$sample_dt) >= 2008)

df.sf.NWIS.keep<-df.sf.NWIS.keep%>%filter(Name %in% unique(temp$site_no))

# remove the Long island sites:

df.sf.NWIS.keep<-df.sf.NWIS.keep%>%filter(!Name %in% c("01304000", "01305000", "01304500"))

# map of this new set of sites:

# mapview(df.sf.NWIS.keep)

# filter the predictor set to these sites:

df.NWIS.Predictors<-filter(df.NWIS.Predictors, Name %in% df.sf.NWIS.keep$Name)

##### Correlations

# in the first pass, all predictors were included
# in a second pass, predictors with lots of zeros are removed:

# second pass filter:
# select only numeric columns in df.NWIS.predictors:
x<-df.NWIS.Predictors%>%dplyr::select(where(is.numeric))
#set NAs to zero:
x[is.na(x)] <- 0
# determine the percentage of cells which are zero in eachcolumn, 
# turn into a dataframe, 
# pivot longer,
# filter for predictors with greater than 50% of observations being zero:
x<-lapply(x, function(x){ length(which(x==0))/length(x)})%>%
  bind_rows(.)%>%
  pivot_longer(cols = everything(), values_to = 'Value', names_to = 'Type')%>%
  filter(Value > .5)

# set up a datatframe to feed into eachOLSand sens pipes:

# combined the predictors df with the cq data:

temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.sf.NWIS.keep$Name)%>%
  left_join(., df.NWIS.Predictors, by = 'Name')%>%
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
# for a second pass, deselect the columns names in x:

df.OLS_Sens<-left_join(df.OLS_Sens, df.NWIS.Predictors%>%select(-x$Type), by = 'Name')

# now run correlations between intercepts and slopes and watershed characteristics. to do this: (I orginally did this workflow using n_months (C:\PhD\Research\Mohawk\Code\Mohawk_Regression-analyizing_predictor_df.R)

# set up variablefor number of sites:

n_sites<-dim(df.OLS_Sens)[1] 

# use the corrr package to correlate() and focus() on your variable of choice

df.cor <- df.OLS_Sens %>% 
  correlate() %>% 
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S))%>%
  pivot_longer(cols= c(2:5), names_to = 'CQ_Parameter', values_to = 'Pearson_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Pearson_Correlation*sqrt((n_sites-2)/(1-(Pearson_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero

# then plotresults: todo this:

# create a list of each CQ parameter (4: OLS and Sens slope and intercept)and format it for ggplotting:

l.cor<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>%mutate(term = factor(term, levels = unique(term[order(Pearson_Correlation)])))%>%filter(!between(Pearson_Correlation, -0.25,.25))) # Order by correlation strength

plist<-lapply(l.cor, \(i) i%>%ggplot(aes(x = term, y = Pearson_Correlation, color = sig_0.05)) +
                geom_bar(stat = "identity") +
                facet_wrap('CQ_Parameter')+
                ylab(paste('Pearson Correlation')) +
                xlab("Watershed Attribute")+
                theme(axis.text.x=element_text(angle=40,hjust=1))+
                theme(legend.position="bottom"))

ggpubr::ggarrange(plotlist = plist, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

# univariate plots of the top correlates/small land use types:

# Lets do plots of OLS slope: to do this:

# determine the top correlates:

OLS<-l.cor[[2]]%>%arrange(desc(Pearson_Correlation))

# make univariate plots (facets) of different land uses and OLS slopes:

df.OLS%>%left_join(.,df.NWIS.Predictors%>%select(c(Name, OLS$term)), by = 'Name')%>%
  pivot_longer(cols = 4:10, names_to = 'Type', values_to = 'Value')%>%
  drop_na(Value)%>%
  ggplot(., aes(x = Value, y = OLS.S))+
  geom_smooth(method = 'lm')+
  geom_point()+
  facet_wrap('Type', scales = 'free')


# CQplot of all sites:

ggplot(temp, aes(x = log_Q, y = log_C))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap('n_sample_rank')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
  
# plot that shows slope and intercept on a map so any spatial correlation would show up?

# merge the lat longs for the sites to df.OLS:

df.OLS%>%
  left_join(., df.NWIS.Predictors%>%select(Name, geometry), by = 'Name')%>%
  st_as_sf(.)%>%
  mapview(., zcol = 'OLS.S')

#### Grouping CQ curves (stationary, mobilization, dilutionary, complex) ####

# create a list of dataframes for each sites CQ observations:

temp<-df.NWIS.TP_CQ%>%
  rename(Name = site_no)%>%
  filter(Name %in% df.sf.NWIS.keep$Name)%>%
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

ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
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

ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
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

# based on this plot, I would chose the following sites as complex:

complex_sites<-unique(df_Seg.2$site)[c(3,18,24,25,26,34,35,42,43)]

df_Seg.2<-mutate(df_Seg.2, Type = ifelse(site %in% complex_sites, 'Complex', Type))

# make plot:

ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("purple", "red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  geom_line(aes(x = Q, y = Seg_C), color = 'yellow', size = 1.5)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# add the type to the mapping df:

df.sf.NWIS.keep.2<-left_join(df.sf.NWIS.keep, distinct(df_Seg.2, site, .keep_all = T)%>%select(site, Type, n_sample_rank), by = c('Name'='site'))%>%
  select(Name, Type, n_sample_rank)%>%
  arrange(n_sample_rank)%>%
  mutate(n_sample_rank=1:nrow(.))%>%
  left_join(.,df.NWIS.TP%>%group_by(site_no)%>%summarise(n=n()), by = c('Name'='site_no'))%>%
  mutate(NEW = case_when(n < 50 ~ .05,
                         n >=50 & n < 100 ~ .25,
                         n >=100 & n < 500 ~ .5,
                         n >=500 ~ .9))
# map:

# mapview(df.sf.NWIS.keep.2, zcol = 'Type', alpha.regions = 'NEW')

#### Categorizing land use ####

# USGS criteria:
# Agricultural sites have >50% agricultural land and ≤5% urban land;
# urban sites have >25% urban and ≤25% agricultural land; 
# undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land; 
# all other combinations of urban, agricultural, and undeveloped lands are classified as mixed

# merge land use from df.NWIS.CDL to df.sf.NWIS.keep.2:

df.sf.NWIS.keep.2<-left_join(df.sf.NWIS.keep.2, df.NWIS.CDL, by = 'Name')

# create the land use class column based on USGS critiera:

df.sf.NWIS.keep.2<-df.sf.NWIS.keep.2%>%
  mutate(USGS.LU = 'Mixed')%>%
  mutate(USGS.LU = case_when(.default = 'Mixed',
    Ag > .50 & Developed <= .05 ~ 'Agriculture',
    Developed > .25 & Ag <= .25 ~ 'Urban',
    Developed <= .05 & Ag <= .25 ~ 'Undeveloped')
  )

# merge the OLS slopes with this df:

df.sf.NWIS.keep.2<-left_join(df.sf.NWIS.keep.2, df.OLS_Sens[,1:5], by = 'Name')

# pivot longer for geom_box + facet:

data<-df.sf.NWIS.keep.2%>%
  pivot_longer(cols = 15:18, names_to = 'CQ_parameter', values_to = 'Value')%>%
  mutate(USGS.LU=factor(USGS.LU))

# make ggplot:

my_xlab <- paste(levels(factor(df.sf.NWIS.keep.2$USGS.LU)),"\n(N=",table(factor(df.sf.NWIS.keep.2$USGS.LU)),")",sep="")

ggplot(data, aes(x=USGS.LU, y=Value))+
  geom_boxplot(varwidth = TRUE, alpha=0.2)+
  scale_x_discrete(labels=my_xlab)+
  facet_wrap('CQ_parameter', scales = 'free')+
  stat_compare_means(method = "anova", label.y = 2)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0.5")  
  

# finally save image to workspace:

# save.image(file = 'Processed_Data/NWIS.Rdata')

load('Processed_Data/NWIS.Rdata')




