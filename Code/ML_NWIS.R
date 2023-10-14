# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

library(climateR)
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

####################### Load in Prior data #######################

load('Processed_Data/NWIS.Rdata')

# keep only the datasets needed for this script:

rm(list=setdiff(ls(), c("df.NWIS.TP_CQ", "df.G2", "df.sf.NWIS.keep", "df.NWIS.Q")))

####################### Functions #######################

source("Code/Ryan_functions.R")

####################### Goal of code #######################

# 

####################### Workflow #######################

# come up with a dataset:
# Gauges 2 sites in NYS that have TP data that have a time series of TP data:

####----Diagnostics to filter sites----####

# filter the raw TpCQ df to the 89 gauges 2 sites with >20 TP CQ observations

d1<-df.NWIS.TP_CQ%>%filter(site_no %in% df.G2$STAID)

# add a column for the number of years of C data: to do this:

# create a summarized df and filter to the sites with over 10 years of data:

temp<-d1%>%
  group_by(site_no)%>%
  summarize(n_years = max(year(sample_dt))-min(year(sample_dt)))%>%
  filter(n_years >= 10)

# filter d1 to just these sites:

d1<-d1%>%filter(site_no %in% temp$site_no)

# plot the time series of these sites:

# ggplot(d1, aes(x = sample_dt, y = result_va))+
#   geom_point()+
#   facet_wrap('n_sample_rank')+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )

# I want to look at a map of the coefficient of variation of the TP samples:
  
# calculate the CV of samples:

TP<-d1%>%group_by(site_no)%>%summarize(cv = cv(result_va, na.rm = T))

# merge the watershed boundaries to this df (need to also convert to sf)

TP<-left_join(TP, df.sf.NWIS.keep%>%select(Name, geometry), by = c('site_no'='Name'))%>%
  st_as_sf(.)%>%
  arrange(desc(cv))%>%
  mutate(cv_rank = 1:nrow(.))

# map:

# mapview(TP, zcol = 'cv')

# I want to look at the time series plots ordered by cv:

# add the cv column to the rawtpdata df:

d1<-left_join(d1, TP%>%st_set_geometry(., NULL), by = 'site_no')

# ggplot(d1, aes(x = sample_dt, y = result_va))+
#   geom_point()+
#   facet_wrap('cv_rank')+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )

# diagnostic the cv:

x<-d1[d1$site_no=="01422747",]

# ggplot(x, aes(x = sample_dt, y = result_va))+
#   geom_point()

x<-d1[d1$site_no=="0423205010",]

# ggplot(x, aes(x = sample_dt, y = result_va))+
#   geom_point()

####----Feature Selection----####

# in the first pass, spearman rank correlaitons showed that there was 
# high correlaitons with the annaul precip values in the database. So I removed those as well as 
# other annually tabulated values I demeed necessary.
# the seocnd pass uses the edited gauges 2 predictorset:

# read in the edited G2 dataset and convert to df:

df.G2_edited <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm_EDITED.xlsx")

df.G2_edited[[24]]<-NULL

df.G2_edited<-reduce(df.G2_edited, full_join, by = "STAID")

# merge the raw TP CQ data with its gauges 2 info:

d1<-left_join(d1, df.G2_edited, by = c('site_no'='STAID'))

# perform spearman rank correlaitons: to do this:

# filter to just numeric columns:

temp<-d1%>%select_if(is.numeric)

# remove columns with all NAs:

temp<-temp[colSums(!is.na(temp)) > 0]

# perform cor of dataframe:

temp<-as.data.frame(cor(temp[-1], temp[1],use="complete.obs", method = 'spearman'))

# create column of row names (which are thepredictor names) and set the row names to 1:n:

temp$Predictor<-rownames(temp)

rownames(temp)<-1:nrow(temp)

# arrange by abs(spearman rank correlation):

temp<-arrange(temp, desc(abs(result_va)))

# create a seperate df for the top 35 and add a rankcolumn:

temp<-temp[1:35,]%>%
  arrange(desc(result_va))%>%
  mutate(Spear_Rank = 1:35)%>%
  mutate(Spear_Rank_name = paste(Spear_Rank, Predictor))

# univariate plots:

# d1%>%
#   select(c(site_no, result_va, temp$Predictor))%>% # select the top 35 predictors (in temp)
#   pivot_longer(cols = -c(site_no, result_va), names_to = 'Predictor', values_to = 'Value')%>% # pivot longer to aid ggplot
#   mutate(across(Predictor, factor, levels=temp$Predictor))%>% # add the factor levels based on a rank esabilished in temp to get facets to plot in descending spearman rank value order. 
#   ggplot(., aes(x = Value, y = result_va))+
#   geom_point()+
#   geom_smooth(method = 'lm')+
#   facet_wrap('Predictor', scales = 'free')

####----Adding Features----####

# I want to add the time sensitive predictors to thedataset:
# these are antcedent flow, rainfall, and temperature (I will use daily mean temps for all temp variables) 
# I like the lags and deltas used in Harrison et al (read it in FOR 797)
# and I like the cummalitve used in Adedeji et al.
# I think rainfall will capture most of the concerns with flow?
# Im just going to pick a few and run with it:

# Flows: 

# flow on day of, and 1,2,3,4,5,10,14,21 days prior to the sample
# delta flow between day of sample and 1,2,3,4,5, 10, 14, 21, days prior to sample

# Daily flows are already apired with the CQ observations
# To do antecdent flows the raw daily data is needed. This as already downloaded for the TP sites in NWIS.R

# filter all TP sites rawdaily flow data to the 55 sites in d1:

df.ML_NWIS.Q<-df.NWIS.Q%>%filter(site_no %in% d1$site_no)

# Each site+date is going to have 8+8 (lag+delta) = 16 flow features

# initalize a dataframe:

lag_deltas<-c(0,1,2,3,4,5,10,14,21)

lag_deltas_char<-as.character(lag_deltas)

lag_names<-paste0('Flow_lag_', lag_deltas_char, '_day')

delta_names<-paste0('Flow_delta_', lag_deltas_char, '_day')

columnnames<-c('site_no', 'Date', lag_names, delta_names)

df.ML_NWIS.Q_features<-setNames(data.frame(matrix(ncol = length(x), nrow = 1)), columnnames)

# loop through the sites and extract the features:

for(i in seq_along(unique(d1$site_no))){
  
  # filter the flows for a site:
  
  df.Q<-filter(df.ML_NWIS.Q, site_no == unique(d1$site_no)[i])
  
  # filter the samples for the site:
  
  df.C<-filter(d1, site_no == unique(d1$site_no)[i])
  
  # filter this down to unique dates. If there are more than one sample on a given date, take the mean:
  
  df.C<-df.C%>%group_by(sample_dt)%>%summarise(result_va_mean = mean(result_va))
  
  # loop through the samples dates:
  
  for (j in seq_along(df.C$sample_dt)){
    
    # determine 0,1,2,3,...21 days prior to sample date:
    
    dates<-as.data.frame(df.C$sample_dt[j]-lag_deltas)%>%rename(Date = 1)
    
    # extract subset of flow dataframe for the dates upto the 21 days prior:
    
    Q.j <- left_join(dates, df.Q, by = 'Date')
    
    # calcualte deltas:
    
    delta.j<-c(0,Q.j$X_00060_00003[-1]-Q.j$X_00060_00003[1])
    
    # Make wide dataframes of the lags and delta
      
    lag.j<-data.frame(site_no = unique(d1$site_no)[i], Date = df.C$sample_dt[j], lag = Q.j$X_00060_00003, name = lag_names) %>%
      pivot_wider(names_from = name, values_from = lag)
    
    delta.j<-data.frame(site_no = unique(d1$site_no)[i], Date = df.C$sample_dt[j], delta = delta.j, name = delta_names) %>%
      pivot_wider(names_from = name, values_from = delta)
    
    # merge these two:
    
    df.j<-left_join(lag.j,delta.j,by = c('site_no', 'Date'))
    
    # append to outside df:
    
    df.ML_NWIS.Q_features<-bind_rows(df.ML_NWIS.Q_features, df.j)
    
  }
  
}

# checkthis


# Rainfall:

# 1(24 hr),2,3,4,5,10,14 and 21 day cumulative rainfall prior to sample date







save.image(file = 'Processed_Data/ML_NWIS.Rdata')
  
# load("Processed_Data/ML_NWIS.Rdata")






