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

ncode<-'TN'






















#### NWIS Query ####

# this code adapted from NYS_site_ranking.R

# Download the metadata for sites with daily flow data. To do this:
# alread did this in Code/NWIS.R:

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
  mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))%>%
  filter(nflowdays>=20)

# NOTE: whatNWISsites(stateCd = "NY", parameterCd = "00060") returns a dataframe with rnow = 1465
# while whatNWISdata(stateCd = "NY", parameterCd = "00060") returns a dataframe with nrow = 2249 ?!?!
# this is because there are duplicates. Using the group_by, slice functions gives a dataframe with
# nrow = 1452 for the whatNWISdata thing.

# From the sites with daily flow data, determine which sites also have over 99 discrete samples . To do this:

# use function I created (sourced) to get a dataframe of sites for just TN with #samples = 20 threshold:

# df.NWIS.TN_sites<-fun.df.Pair_consit_flow(pcode, df.NWIS.Q_sites, n_samples = 20, state = 'NY')

# write.csv(df.NWIS.TN_sites, "Raw_Data/df.NWIS.TN_sites.csv", row.names=FALSE)

df.NWIS.TN_sites<-read.csv("Raw_Data/df.NWIS.TN_sites.csv", colClasses = c(site_no = "character"))

# download the raw daily flow data for these sites. To do this:
# I already downloaded the daily flow data for the 102 TN sites
# which ever sites do notoverlap with TN will be downloaded here. To do this:

# read in the flow data for the TN sites:

# df.NWIS.Q<-read.csv("Raw_Data/df.NWIS.Q.csv", colClasses = c(site_no = "character"))%>%mutate(Date = as.Date(Date))

# determine which of the TN sites dont overlap:

# df.overlap<-df.NWIS.TN_sites%>%filter(!site_no %in% df.NWIS.Q$site_no)

# download the raw flow data for the sites still needed:

# df.NWIS.Q.TN_needed<-readNWISdv(siteNumbers = df.overlap$site_no, parameterCd = '00060', startDate = "", endDate = "", statCd = "00003")

# unique(df.NWIS.Q.TN_needed$site_no) # only get three of 25 that need... sucks but idk

# merge with the TN sites flow data:

# df.NWIS.Q.combined<-bind_rows(df.NWIS.Q,df.NWIS.Q.TN_needed)

# write out this df to file:
  
# write.csv(df.NWIS.Q.combined, "Raw_Data/df.NWIS.Q.combined.csv", row.names=FALSE)

# read it back in as 'df.NWIS.Q.for_TN':

df.NWIS.Q.for_TN<-read.csv("Raw_Data/df.NWIS.Q.combined.csv", colClasses = c(site_no = "character"))%>%mutate(Date = as.Date(Date))

# filter to just the TN sites:

df.NWIS.Q.for_TN<-df.NWIS.Q.for_TN%>%filter(site_no %in% df.NWIS.TN_sites$site_no)

# download the raw discrete TN sample data:

# df.NWIS.TN<-readNWISqw(siteNumbers = df.NWIS.TN_sites$site_no, parameterCd = pcode)

# write.csv(df.NWIS.TN, "Raw_Data/df.NWIS.TN.csv", row.names=FALSE)

df.NWIS.TN<-read.csv("Raw_Data/df.NWIS.TN.csv", colClasses = c(site_no = "character"))%>%mutate(sample_dt = as.Date(sample_dt))

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

# df.NWIS.TN_site_metadata<-readNWISsite(siteNumbers = unique(df.NWIS.TN_CQ$site_no))

# write.csv(df.NWIS.TN_site_metadata, "Raw_Data/df.NWIS.TN_site_metadata.csv", row.names=FALSE)

df.NWIS.TN_site_metadata<-read.csv("Raw_Data/df.NWIS.TN_site_metadata.csv", colClasses = c(site_no = "character"))

# then select just the site number and DA column:

df.DA<-df.NWIS.TN_site_metadata%>%
  select(site_no, drain_area_va)

# finally merge with df.NWIS.TN_CQ and create a new Q column with area normalized flows (not worrying about units right now): 
# Note: will filter for NA in C and Q for breakpoint analysis in the next step as to keep the full list of sites with CQ pairs in this dataframe.
# Note: Some sites returned NA on draiange areas in readNWISsite, but I'll delinate anyways so I want the full list:

df.NWIS.TN_CQ<-left_join(df.NWIS.TN_CQ, df.DA, by = 'site_no')%>%
  mutate(Q_yield = X_00060_00003/drain_area_va)

#









































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

























#### Apply filters to sites ####

# the first filter is to remove data points and thus potentially entire sites
# that are prior to 2000:

temp1<-df.NWIS.TN_CQ%>%filter(sample_dt >= as.Date('2001-01-10'))

unique(temp1$site_no)

# 39 sites, BUT:

# how many of these sites have 20 paired CQ observations?

temp1.1<-temp1%>%
  rename(Name = site_no)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))%>%
  group_by(Name)%>%
  summarise(n=n())%>%
  filter(n>=20)

dim(temp1.1)[1] 

# only 31

# reduce temp1 to these sites:

temp1<-temp1%>%filter(site_no %in% temp1.1$Name)

# the next filter is to remove sites that are on long island
# to do this use latitude:

temp2<-filter(df.NWIS.TN_site_metadata, dec_lat_va >40.9364)

# now use this site list to filter down the df.TN_CQ:

temp3<-temp1%>%filter(site_no %in% temp2$site_no)

unique(temp3$site_no)

# 31 sites

# the last filter is gauges 2:

# read in all sheets using function

# l.G2 <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm.xlsx")

# remove the last element (does notcomtain useful info)

# l.G2[[27]]<-NULL

# and convert to a single df

# df.G2<-reduce(l.G2, full_join, by = "STAID")

# save df.G2:

# save(df.G2, l.G2, file = 'Processed_Data/df.G2.Rdata')

# load df.G2:

load('Processed_Data/df.G2.Rdata')

# filter on G2:

temp4<-temp3%>%filter(site_no %in% df.G2$STAID)

unique(temp4$site_no)

# 17 sites

# save this as the df to move onto future analysis:
# need to alsoremove samples after 2001:

df.TN_CQ<-temp4%>%filter(sample_dt >= as.Date('2001-01-10'))

# map:

# df.NWIS.TN_site_metadata%>%
#   filter(site_no %in% temp4$site_no)%>%
#   st_as_sf(.,coords=c('dec_long_va','dec_lat_va'), crs = 4326, remove = FALSE)%>%
#   mapview(.)

# number of sites in gauges 2:

x<-filter(df.NWIS.TN_site_metadata, site_no %in%df.G2$STAID)

# 68 sites

# number of sites in gauges 2 and not on LI (i.e. without post 2001 filter):

x<-temp2%>%filter(site_no%in% df.G2$STAID)

# 53 sites























#### Export data ####

# save(df.TN_CQ,file = 'Processed_Data/TN.Rdata')























#### Fitting Breakpoints to CQ curves #### 

# create empty list to hold the results of the segmented function (which does the breakpoint analysis)

l_Seg<-list()

# create a matrix to hold the results of the davies test, which determines if a two slope model is warrented over a single slope model:

davies.test.matrix<-NULL

# create a new dataframe of only paired CQ observations (such that the breakpoit analysis function runs smoothly) (I didnt want to do this in loop for some reason, I cant remeber why but it wouldnt work):

df.TN_CQ_for_BP<-df.TN_CQ%>%
  drop_na(result_va, Q_yield)

# create a vector of ordered unique site names:

temp.loop<-sort(unique(df.TN_CQ_for_BP$site_no))

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
    
    df<-df.TN_CQ_for_BP%>%
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

# p<-ggplot(df_Seg, aes(x = log(Q_real), y = log(C)))+
#   geom_point()+
#   geom_smooth(method = 'lm')+
#   geom_line(aes(x = Q, y = Seg_C), color = 'tomato')+
#   facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )

# p


# one last thing: lets look at a map of these:

# map.NWIS.TN_sites<-df.NWIS.TN_site_metadata%>%
#   rename(longitude=8,latitude=7)%>%
#   drop_na(latitude,longitude)%>%
#   st_as_sf(.,coords=c('longitude','latitude'), crs = 4326)%>%
#   mutate(site_id= paste(site_no, station_nm), n = NA, .before = 1)%>%
#   select(c(1:3))

# mapview(map.NWIS.TN_sites)




























#### Calculating Average Annual Consituent Yield ####
# (for use in correlaitons and MLR)

# step 1. calcualte average annual hydrograph for each site. to do this:

# filter the daily flow records to just the sites used:

l.avg_ann_hydro<-df.NWIS.Q.for_TN%>%filter(site_no %in% df.TN_CQ$site_no)

# split the daily flow records into list of dfs

l.avg_ann_hydro<-split(l.avg_ann_hydro, f = l.avg_ann_hydro$site_no)

# strip out year from date (just converts it to 2023) for each dataframe, then group by and summarize to get mean annual flow for each day of calendar year:

l.avg_ann_hydro<-lapply(l.avg_ann_hydro, \(i) i%>%mutate(Date = as.Date(Date))%>%mutate(year = year(Date), Date = as.Date(format(Date, "%m-%d"), "%m-%d"))%>%group_by(Date)%>%summarize(mean_Q = mean(X_00060_00003, na.rm= T)))

# step 2. predict C for each of day of the year using the mean flow rate for each site in the list:todothis:

# create a list of linear model objects for each sites CQ relaitonship:

l.C_daily<-df.TN_CQ%>%
  rename(Name = site_no)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))%>%
  split(., f =.$Name)%>%
  lapply(., \(i) lm(log_C ~ log_Q, i))

# make sure the order of the sites in the flow and lm lists are the same:

names(l.C_daily)==names(l.avg_ann_hydro)

# use this list to create list of bias correction factors:

l.BCF<-lapply(l.C_daily, \(i) exp((summary(i)$sigma^2)/2))

# predict the concentration for each day for each site: this predicts the concnetration in log space, so transforming into real space here as well using BCF:

l.C_daily<-lapply(names(l.C_daily), \(i) l.BCF[[i]]*exp(as.numeric(predict(l.C_daily[[i]], data.frame(log_Q = log(l.avg_ann_hydro[[i]]$mean_Q) ) ) ) ) )%>%purrr::set_names(names(l.avg_ann_hydro))

# step 3. multiply predicted average daily C by daily Q to get daily load, sum up, then convert to yield. To do this:

# save just the daily Q:

l.Q.temp<-lapply(l.avg_ann_hydro, \(i) i$mean_Q)

# set up list of drainage areas:

l.DA<-readNWISsite(siteNumbers = names(l.C_daily))%>%select(site_no, drain_area_va)%>%split(., .$site_no)%>%lapply(., \(i) i$drain_area_va)

# multiply daily C by daily Q and convert units:
# mg/L * ft^3/sec * 28.3168 L/ft3 * 86400 sec /day * 1g/1000mg * 1kg/1000g * 1/mi^2 * 1mi^2/258.999 ha = kg/ha/day

unit_conversion=28.3168*86400*(1/1000)*(1/1000)*(1/258.999)

l.Yield<-lapply(names(l.C_daily), \(i) l.C_daily[[i]]*l.Q.temp[[i]]*(1/l.DA[[i]])*unit_conversion)%>%purrr::set_names(names(l.avg_ann_hydro))

# sum up for the year and convert to df

df.Yield <- lapply(l.Yield, sum)%>%dplyr::bind_rows(., .id = 'Name')%>%pivot_longer(cols = everything(), names_to = 'Name', values_to = 'Yield') # units of kg/ha/year

#























#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Correlations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# reduce the predictor set down to those we believe are defensible:

# I also went in to the variable description excel sheet and isolated these variables as well:

temp1<-read.csv('Processed_Data/G2_pred_to_keep_2.csv')[,1]

# combining these with the vairbales in entire sheets as well:

pred_to_keep<-c("HYDRO_DISTURB_INDX", 
                "pre1990_STOR", 
                names(l.G2$Landscape_Pat), 
                names(l.G2$LC06_Basin),
                names(l.G2$LC_Crops),
                temp1,
                'HGA', 'HGB', 'HGC', 'HGD',
                "ELEV_MEDIAN_M_BASIN",
                "ELEV_STD_M_BASIN",
                "RRMEDIAN"
)

# removing 'STADID':

pred_to_keep <- pred_to_keep[!(pred_to_keep %in% "STAID")]

# removing the duplicates of DEV and PDEN

pred_to_keep <- pred_to_keep[!(pred_to_keep %in% c("DEVOPENNLCD06","DEVLOWNLCD06","DEVMEDNLCD06","DEVHINLCD06", "PDEN_DAY_LANDSCAN_2007", "PDEN_NIGHT_LANDSCAN_2007"))]

# filtering df.G2 to this predictor list:

df.G2.reduced<-df.G2%>%select(c('STAID', pred_to_keep))

# take the average of the RIP and MAIN for each land use:

df.G2.reduced<-df.G2.reduced%>%filter(STAID %in% df.TN_CQ$site_no)%>%
  rowwise() %>%
  mutate(MAIN_RIP_DEV_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("DEV")), na.rm = TRUE),
         MAIN_RIP_FOREST_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("FOREST")), na.rm = TRUE),
         MAIN_RIP_PLANT_avg = mean(c_across(starts_with("MAINS") | starts_with("RIP") & ends_with("PLANT")), na.rm = TRUE),
         .keep = 'unused'
  )

# adding up the other crops:

df.G2.reduced<-df.G2.reduced%>%
  mutate(sum_of_misc_crops = sum(c_across(ends_with(c("COTTON", "RICE", "SORGHUM", "SUNFLOWERS", "PEANUTS", "BARLEY", "DURUM_WHEAT", "OATS", "DRY_BEANS", "POTATOES", "ORANGES", 'OTHER_CROPS', 'IDLE'))), na.rm = TRUE), .keep = 'unused')

names(df.G2.reduced)

# combined the reduced predictors df with the already filtered CQ data:

df.Correlation<-df.TN_CQ%>%
  rename(Name = site_no)%>%
  select(Name, sample_dt,result_va, X_00060_00003)%>%mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))%>%
  left_join(., df.G2.reduced, by = c('Name'='STAID'))

# create a dataframe with OLS and sens slope intercept and slopes:

df.OLS<-df.Correlation%>%
  group_by(Name)%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.I = OLS.co[1], 
            OLS.S = OLS.co[2])
  }) %>%
  ungroup

df.Sens<-df.Correlation%>%
  group_by(Name)%>%
  do({ Sens.co<-zyp.sen(log_C~log_Q,.)
  summarize(., Sen.I = Sens.co$coefficients[[1]],
            Sen.S= Sens.co$coefficients[[2]])
  }) %>%
  ungroup

# merge OLS and Sens:

df.OLS_Sens<-left_join(df.OLS, df.Sens, by = 'Name')

# merge constiuent yield:

df.OLS_Sens<-left_join(df.OLS_Sens, df.Yield, by = 'Name')

# merge back the watershed characteristic data to this dataframe:

df.OLS_Sens<-left_join(df.OLS_Sens, df.G2.reduced, by = c('Name'='STAID'))

# now run correlations between intercepts and slopes and watershed characteristics. to do this: (I orginally did this workflow using n_months (C:\PhD\Research\Mohawk\Code\Mohawk_Regression-analyizing_predictor_df.R)

# set up variable for number of sites:

n_sites<-dim(df.OLS_Sens)[1] 

# use the corrr package to correlate() and focus() on your variable of choice.

df.cor <- df.OLS_Sens %>% 
  correlate(method = 'spearman') %>%
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S, Yield))%>%
  pivot_longer(cols= c(2:6), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero

# I also want to try standardizing the predictors:

df.cor.0to1<-df.OLS_Sens %>% 
  mutate(across(7:last_col(), ~ (.-min(.))/(max(.)-min(.))))%>% # could add term: .names = "{.col}_standarized")
  correlate(method = 'spearman') %>%
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S, Yield))%>%
  pivot_longer(cols= c(2:6), names_to = 'CQ_Parameter', values_to = 'Spearman_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Spearman_Correlation*sqrt((n_sites-2)/(1-(Spearman_Correlation)^2))), n_sites-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero

# looking at the regualr and 0 to 1 standarized spearman correlaitons, they are the same

df.cor$Spearman_Correlation==df.cor.0to1$Spearman_Correlation

# so not going to proceed with 0 to 1

# then plot results: to do this:

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
                scale_color_manual(values = c("not" = "red", "sig" = "blue"),na.value = NA, drop = FALSE)+
                facet_wrap('CQ_Parameter')+
                ylab(paste('Spearman Correlation')) +
                xlab("Watershed Attribute")+
                theme(axis.text.x=element_text(angle=40,hjust=1))+
                theme(legend.position="bottom"))

# arrange plots:

p1<-ggpubr::ggarrange(plotlist = plist, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

# add title:

p1<-annotate_figure(p1, top = text_grob(paste(ncode, "Correlated Against Gauges 2"),
                                        color = "red", face = "bold", size = 14))
# plot:

p1

# save(p1, file = 'Processed_Data/p1.Rdata')

# univariate plots of the top correlates/small land use types:

# Lets do plots of OLS intercept: to do this:

# determine the top correlates (the numberof the list determined which CQ parameter)
# here 1 = OLS intercept:

OLS<-l.cor[[1]]%>%arrange(desc(Spearman_Correlation))

# make univariate plots (facets) of Y (determined in [[]]) above and predictors: 
# note the predictor column isturned into an ordered factor based on thespearman correlation values
# to makethe facet plots in order:

df.OLS_Sens%>%
  select(Name, OLS$CQ_Parameter[1], OLS$term)%>%
  pivot_longer(cols = c(3:last_col()), names_to = 'Type', values_to = 'Value')%>%
  drop_na(Value)%>%
  # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
  left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
  mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
  ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
  geom_smooth(method = 'lm')+
  geom_point()+
  facet_wrap('Type', scales = 'free')+
  ggtitle(paste('Univariate plots of', OLS$CQ_Parameter[1], 'against Top Gauges 2 Correlates for', OLS$CQ_Parameter[1], 'ordered from highest to lowest Spearman Rank Correlation' ))

# create a function to use lapply to plot all 4 univariate sets:

# number<-5

make_plot<-function(number){
  
  OLS<-l.cor[[number]]%>%arrange(desc(Spearman_Correlation))
  
  x<-df.OLS_Sens%>%
    select(Name, OLS$CQ_Parameter[1], OLS$term)%>%
    pivot_longer(cols = c(3:last_col()), names_to = 'Type', values_to = 'Value')%>%
    drop_na(Value)%>%
    # mutate_if(is.numeric, ~replace(., . == 0, NA))%>%
    left_join(., OLS%>%select(term, Spearman_Correlation), by = c('Type'='term'))%>%
    mutate(Type = factor(Type, levels=unique(Type[order(-Spearman_Correlation,Type)]), ordered=TRUE))%>%
    ggplot(., aes(x = Value, y = !!sym(OLS$CQ_Parameter[1])))+
    geom_smooth(method = 'lm')+
    geom_point()+
    facet_wrap('Type', scales = 'free')+
    ggtitle(paste('Univariate plots of', OLS$CQ_Parameter[1], 'against Top Gauges 2 Correlates for', OLS$CQ_Parameter[1], 'ordered from highest to lowest Spearman Rank Correlation' ))
  
  x
  
}

# use function in lapply (clear plot list first):

# lapply(c(1:5), \(i) make_plot(i))

#






















#### Categorizing land use ####

# USGS criteria:
# Agricultural sites have >50% agricultural land and ≤5% urban land;
# urban sites have >25% urban and ≤25% agricultural land; 
# undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land; 
# all other combinations of urban, agricultural, and undeveloped lands are classified as mixed

# initate a dataframe with the landuse values from gauges 2:

df.NLCD06<-distinct(df.OLS_Sens, Name, .keep_all = T)

# I confirmed that pasture + crops = plant:

(df.NLCD06$PASTURENLCD06+df.NLCD06$CROPSNLCD06)==df.NLCD06$PLANTNLCD06

# and that the ones that start with DEV sum up to the first DEV one.

# create the land use class column based on USGS critiera:
# adjusted thresholds:

df.NLCD06<-df.NLCD06%>%
  mutate(USGS.LU.Adjusted = case_when(.default = 'Mixed',
                                      PLANTNLCD06 > 30 & DEVNLCD06 <= 10 ~ 'Agriculture',
                                      DEVNLCD06 > 10 & PLANTNLCD06 <= 30 ~ 'Urban',
                                      DEVNLCD06+PLANTNLCD06 <= 10 ~ 'Undeveloped'))

# frequency table:

table(df.NLCD06$USGS.LU.Adjusted)

# boxplots of

p<-df.NLCD06%>%
  pivot_longer(cols = 2:5, names_to = 'CQ_parameter', values_to = 'Value')%>%
  ggplot(., aes(x=USGS.LU.Adjusted, y=Value, color =USGS.LU.Adjusted ))+
  geom_boxplot(varwidth = TRUE, alpha=0.2)+
  # scale_x_discrete(labels=my_xlab)+
  facet_wrap('CQ_parameter', scales = 'free')+
  stat_compare_means(method = "anova")+      # Add global p-value # , label.y = max(df.datalayers$Value) ???
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0.5") +
  ggtitle('Adjusted USGS Thresholds using Aggregated Data Layers')

# p

#






































#### Grouping CQ curves (stationary, mobilization, dilutionary, complex) ####

# create a list of dataframes for each sites CQ observations:

df.TN_CQ<-df.TN_CQ%>%
  rename(Name = site_no)%>%
  mutate(log_C = log(result_va), log_Q = log(X_00060_00003), C = result_va, Q = X_00060_00003)%>%
  filter(is.finite(log_C))%>%
  filter(is.finite(log_Q))

l.TN_CQ<-df.TN_CQ%>%
  split(., .$Name)

# create lm models for each site:

l.TN_CQ.lm<-lapply(l.TN_CQ, \(i) lm(log_C~log_Q, data=i))

# save the model coef ad pvals:

TN_CQ.coef<-as.data.frame(bind_rows(lapply(l.TN_CQ.lm, \(i) summary(i)$coefficients[,1]), .id = 'site_no'))%>%
  rename(Intercept = 2, Slope = 3)

# save the pvalues 

TN_CQ.pvals<-as.data.frame(bind_rows(lapply(l.TN_CQ.lm, \(i) summary(i)$coefficients[,4]), .id = 'site_no'))%>%
  rename(Intercept.pval = 2, Slope.pval = 3)

# merge the two dfs:

df.lm<-left_join(TN_CQ.coef,TN_CQ.pvals,by='site_no')

# add column for CQ type:

df.lm<-mutate(df.lm, Type = ifelse(Slope.pval>0.05, 'Stationary', ifelse(Slope>0, 'Mobilization', 'Dilution')))

# merge CQ type labels to the df for plotting

df.TN_CQ<-left_join(df.TN_CQ,df.lm%>%select(site_no, Type),by=c('Name'='site_no'))

# create a df for CQplot of all sites with BP analysis:

df_Seg.2<-filter(df_Seg, site %in% df.lm$site_no)%>%
  left_join(.,df.lm%>%select(site_no, Type),by=c('site'='site_no')) 

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

# p

# looking at this plot I want to add a fourth CQ type for complex, if the slopes of the BP analysis look widely different. 
# I will start with calcuating the angle between pre-post BP slope and see which sites get signled out:

# add a new column with the angle between the two lines:

df_Seg.2<-df_Seg.2%>%
  mutate(slope_angle=factor(round(atan(abs((Slope2-Slope1)/(1+(Slope2*Slope1)))),1)))

# I want to add land use as color as well. Merge the plotting df with the land use df:

df_Seg.2<-left_join(df_Seg.2, df.NLCD06%>%select(Name, USGS.LU.Adjusted), by = c('site'='Name'))

# create color pallete for the slope angle:

hc<-heat.colors(length(unique(df_Seg.2$slope_angle)), rev = T)

# make plot:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
  geom_smooth(method = 'lm')+
  ylab('log(TN)')+
  new_scale_color() +
  geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
  geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
  scale_color_manual(name = "Slope Angle", values = hc)+
  facet_wrap(dplyr::vars(n_sample_rank), scales = 'free')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .35)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","yellow", "green"))


p

# export df_Seg.2:

# save(df_Seg.2, file = 'Processed_Data/TN.df_Seg.2.Rdata')

# # based on this plot, I think I need to zoom in on each one:
# 
# # create a plotting function for each site:
# 
# p.fun<-function(df){
#   ggplot(df, aes(x = log(Q_real), y = log(C)))+
#     geom_point(aes(color = Type))+
#     scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
#     geom_smooth(method = 'lm')+
#     new_scale_color() +
#     geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
#     geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
#     scale_color_manual(name = "Slope Angle", values = hc)+
#     ggtitle(df_Seg.2$site[df_Seg.2$n_sample_rank==i])
# }
# 
# # use function in lappy to make lots of plots (clear plot list first)
# 
# lapply(sort(unique(df_Seg.2$n_sample_rank)), \(i) df_Seg.2%>%filter(n_sample_rank==i)%>%p.fun(.))
# 
# # now looking at the expanded plots for each site, I feel that
# # non really could be justified as complex! it feels like the breakpoint
# # analysis lines aren't 'real'
# 
# # My first thought is to color the CQ points based on season, AMC, etc
# # this is so overwhleming. 
# 
# # I think I'm going to pause here in the analysis, and carry on with the
# # other constituents. Once I have those up to this point I can check in with Chuck and Steve






















#### Steve's Idea: triangle CQ plot ####

# set up dataframe:

# add percent ag and urban to df.lm:

df.tri<-left_join(df.lm, df.NLCD06%>%select(Name, DEVNLCD06, PLANTNLCD06), by = c('site_no'='Name'))%>%select(site_no, Slope, Type, DEVNLCD06, PLANTNLCD06)%>%mutate(Slope = round(Slope, 2), DEVNLCD06=round(DEVNLCD06, 2), PLANTNLCD06=round(PLANTNLCD06, 2))

# create plot:

p<-ggplot(df.tri, aes(x=DEVNLCD06, y=PLANTNLCD06)) +
  geom_point(aes(shape=as.factor(Type), fill=abs(Slope)), size = 4) +
  scale_shape_manual(values = c("Mobilization" = 24, "Dilution" = 25, "Stationary" = 22),
                     guide = guide_legend(override.aes = list(fill = "pink")))+
  scale_fill_gradient(low = "yellow", high = "red")+
  xlab('Percent Developed')+
  ylab('Percent Agriculture')+
  ggtitle(paste(ncode, 'CQ type and slope magnitude as a function of percent Ag and Developed Land'))

p$labels$fill <- "Slope Magnitude"
p$labels$shape <- "CQ Type"

p


#


























#### MLR ####

# create a list of dfs, one for each cQ parameter, with the values of the spearman correlations for each predictor:

l.cor.MLR<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(desc(abs(Spearman_Correlation)))%>%  
           # slice_head(n = 7)%>%
           # bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           # distinct(term, .keep_all = T)%>%
           # filter(!between(Spearman_Correlation, -0.25,.25))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           # mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           # filter(sig_0.05 == 'sig')%>%
           as.data.frame())


x<-l.cor.MLR[[1]]

# create a list of dfs from df.OLS.Sens, which contains the values of the predictors, subseting using the names in each l.cor[[i]]$term:

l.cor.MLR.full<-lapply(1:5, \(i) df.OLS_Sens%>%
                         select(i+1, l.cor.MLR[[i]]$term)%>%
                         as.data.frame()%>%
                         rename(term = 1))

names(l.cor.MLR.full)<-names(df.OLS_Sens)[2:6]

y<-l.cor.MLR.full[[1]]

y.round<-round(y, 2)

# create another list but keep the site names for use at end of loop in site outliers:

l.cor.MLR.full.w_name<-lapply(1:5, \(i) df.OLS_Sens%>%
                                select(1, i+1, l.cor.MLR[[i]]$term)%>%
                                as.data.frame())

names(l.cor.MLR.full.w_name)<-names(df.OLS_Sens)[2:6]

z<-l.cor.MLR.full.w_name[[1]]

# Trying no feature selection:

# loop through the 5 CQ parameters:

i<-1

m.list<-list()

for (i in 1:5){
  
  df<-l.cor.MLR.full[[i]]
  
  library(caret)
  
  set.seed(123)
  
  # Set up repeated k-fold cross-validation
  
  train.control <- trainControl(method = "cv", number = 10)
  
  # Train the model
  
  step.model <- train(term ~., data = df,
                      method = "leapForward", 
                      tuneGrid = data.frame(nvmax = 1:5),
                      trControl = train.control
  )
  
  
  # look at model results:
  
  step.model$results
  step.model$bestTune
  summary(step.model$finalModel)
  coef(step.model$finalModel, as.numeric(step.model$bestTune))
  
  # get df of just the predictors in this model:
  
  pred<-names(coef(step.model$finalModel, as.numeric(step.model$bestTune)))[-1]
  
  df.2<-df%>%select(term, pred)
  
  # make lm:
  
  m<-lm(term~., data=df.2)
  
  m.list[[i]]<-m
  
  # summary(m)
  
  # look at univarite plots:
  
  p<-df.2%>%pivot_longer(cols = 2:last_col(), names_to = 'Predictor', values_to = 'Value')%>%
    rename(!!names(l.cor.MLR.full)[i]:=term)%>%
    ggplot(., aes(x=Value, y=!!sym(names(l.cor.MLR.full)[i])))+
    geom_point()+
    geom_smooth(method = 'lm')+
    facet_wrap('Predictor', scales = 'free')
  
  p
  
  # look at model residuals:
  
  # plot(m)
  
}

# set m.list names:

names(m.list)<-names(l.cor.MLR.full)

tab_model(m.list, dv.labels = names(m.list), title = paste('Comparison of MLR models for',ncode, 'using forward selection implemented in caret::train'), file="temp.html")


#
