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

rm(list=setdiff(ls(), c("df.NWIS.TP_CQ", "df.G2", "df.sf.NWIS", "df.NWIS.Q")))

####################### Functions #######################

source("Code/Ryan_functions.R")

####################### Goal of code #######################

# 

####################### Workflow #######################

# come up with a dataset:
# Gauges 2 sites in NYS that have TP data that have a time series of TP data:

####----Diagnostics to filter sites----####

# how many of the 137 NYS TP CQ >20 sites and how many of the 89 gauges 2 sites delineated correctly:

vect.NWIS<-vect(df.sf.NWIS)
vect.NWIS$area_KM2<-expanse(vect.NWIS, unit="km")
df.sf.NWIS$area_sqmi<-expanse(vect.NWIS, unit="km")*0.386102 # add DA in mi2 to df.sf:
temp<-readNWISsite(df.sf.NWIS$Name) # add NWIS tabulated area to df sf: first download the metadata for the sites using readNWISsite, then merge drain_area_va column to sfdf
df.sf.NWIS<-left_join(df.sf.NWIS, temp[,c(2,30)], by = c('Name'='site_no'))
df.sf.NWIS$Delination_Error<-(df.sf.NWIS$area_sqmi-df.sf.NWIS$drain_area_va)/df.sf.NWIS$drain_area_va # calcuate the percent error of the delination:

df.sf.NWIS.keep<-df.sf.NWIS%>%filter(abs(Delination_Error)<=.02)

# 103

x<-filter(df.sf.NWIS.keep, Name %in% df.G2$STAID)

# 70

# how many of the 137, 89, and 70 have over 10 years of data?

x<-df.NWIS.TP_CQ%>%
  group_by(site_no)%>%
  summarize(n_years = max(year(sample_dt))-min(year(sample_dt)))%>%
  filter(n_years >= 10)

# 66 of 137

x<-df.NWIS.TP_CQ%>%
  filter(site_no %in% df.G2$STAID)%>%
  group_by(site_no)%>%
  summarize(n_years = max(year(sample_dt))-min(year(sample_dt)))%>%
  filter(n_years >= 10)

# 55 of 89

x<-df.NWIS.TP_CQ%>%
  filter(site_no %in% df.G2$STAID)%>%
  filter(site_no %in% df.sf.NWIS.keep$Name)%>%
  group_by(site_no)%>%
  summarize(n_years = max(year(sample_dt))-min(year(sample_dt)))%>%
  filter(n_years >= 10)

# 43 of 70

# map of these 43:

# filter(df.sf.NWIS.keep, Name %in% x$site_no)%>%mapview(.)

filter(df.sf.NWIS.keep, Name %in% x$site_no)%>%st_is_valid()

####----Scenario 1 - the 43 sites----####

# Building 5 feature sets:

# Non-time + time specific features (sites are in gauges 2 and have delineations, i.e. full amount of data)
# non-time specific features + flow (sites are in gauges 2, have flow, but have no delineation)
# Just time specific features (sites are not in gauges 2 but have flow and delinations)
# Just the flow component of the time specific features (sites are not in gauges 2 and have no delineation, but have flow)
# Just non-time specific features (sites are in gauges 2 but have no flow and no delineation - I understand that if a site is in gauges 2 it will have flow, but this situation is mimics when we might recreate the gauges 2 features ourselves)

# I have the non time, so getting the time:

####----Adding Features----####

# I want to add the time sensitive predictors to thedataset:
# these are antcedent flow, rainfall, and temperature (I will use daily mean temps for all temp variables) 
# I like the lags and deltas used in Harrison et al (read it in FOR 797)
# and I like the cummalitve used in Adedeji et al.
# I think rainfall will capture most of the concerns with flow?
# Im just going to pick a few and run with it:

#### Flows: ####

# flow on day of, and 1,2,3,4,5,10,14,21 days prior to the sample
# delta flow between day of sample and 1,2,3,4,5, 10, 14, 21, days prior to sample

# Daily flows are already apired with the CQ observations
# To do antecdent flows the raw daily data is needed. This as already downloaded for the TP sites in NWIS.R

# filter all TP sites rawdaily flow data to the 55 sites in d1:

df.ML_NWIS.Q<-df.NWIS.Q%>%filter(site_no %in% x$site_no)

# Each site+date is going to have 8+8 (lag+delta) = 16 flow features

# initalize a dataframe:

lag_deltas<-c(0,1,2,3,4,5,10,14,21)

lag_deltas_char<-as.character(lag_deltas)

lag_names<-paste0('Flow_lag_', lag_deltas_char, '_day')

delta_names<-paste0('Flow_delta_', lag_deltas_char, '_day')

columnnames<-c('site_no', 'Date', lag_names, delta_names)

df.ML_NWIS.Q_features<-setNames(data.frame(matrix(ncol = length(columnnames), nrow = 1)), columnnames)

# loop through the sites and extract the features:

for(i in seq_along(unique(x$site_no))){
  
  # filter the flows for a site:
  
  df.Q<-filter(df.ML_NWIS.Q, site_no == x$site_no[i])
  
  # filter the samples for the site:
  
  df.C<-filter(df.NWIS.TP_CQ, site_no == x$site_no[i])
  
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
    
    lag.j<-data.frame(site_no = x$site_no[i], Date = df.C$sample_dt[j], lag = Q.j$X_00060_00003, name = lag_names) %>%
      pivot_wider(names_from = name, values_from = lag)
    
    delta.j<-data.frame(site_no = x$site_no[i], Date = df.C$sample_dt[j], delta = delta.j, name = delta_names) %>%
      pivot_wider(names_from = name, values_from = delta)
    
    # merge these two:
    
    df.j<-left_join(lag.j,delta.j,by = c('site_no', 'Date'))
    
    # append to outside df:
    
    df.ML_NWIS.Q_features<-bind_rows(df.ML_NWIS.Q_features, df.j)
    
  }
  
}

# remove first row and flow delta zero day:

df.ML_NWIS.Q_features<-df.ML_NWIS.Q_features[-1,-12]

# save:

# save(df.ML_NWIS.Q_features,file = 'Processed_Data/df.ML_NWIS.Q_features.Rdata')

load('Processed_Data/df.ML_NWIS.Q_features.Rdata')

#### Climate: ####

# 1(24 hr),2,3,4,5,10,14 and 21 day cumulative rainfall prior to sample date
# 0,1(24 hr),2,3,4,5,10,14 and 21 day lag in tmax and tmin:
# 1(24 hr),2,3,4,5,10,14 and 21 day delta in tmax and tmin:

# so the resulting climate df will have:
# 2 (site, date) + 8 (rainfall cummulatives) + 9*2 (tmax and tmin lags) + 8*2 (tmax and tmin deltas)
# = 2+8+9*2+8*2 = 44 columns

# the site+dates are in the previously created df

# create vect of drainage areas:

vect.NWIS<-vect(filter(df.sf.NWIS.keep, Name %in% x$site_no))

# create vectors of climate variables:

vars<-c('pr', 'tmmx', 'tmmn')

# create vectors of functions for group_by and summarize for each climate variables:

vars_funs<-c(`sum`, `max`, `min`)

# initalize a dataframe:

temp_lags<-c(0,1,2,3,4,5,10,14,21); temp_lags_char<-as.character(temp_lags)

rain_cum_and_temp_deltas<-c(1,2,3,4,5,10,14,21); rain_cum_and_temp_deltas_char<-as.character(rain_cum_and_temp_deltas)

rain_colnames<-paste0('Rain_cummulative_', rain_cum_and_temp_deltas_char, '_day')

tmax_lag_colnames<-paste0('Tmax_lag_', temp_lags_char, '_day')

tmin_lag_colnames<-paste0('Tmin_lag_', temp_lags_char, '_day')

tmax_delta_colnames<-paste0('Tmax_delta_', rain_cum_and_temp_deltas, '_day')

tmin_delta_colnames<-paste0('Tmin_delta_', rain_cum_and_temp_deltas, '_day')

columnnames<-c('site_no', 'Date', rain_colnames, tmax_lag_colnames, tmax_delta_colnames, tmin_lag_colnames, tmin_delta_colnames)

df.ML_NWIS.Climate_features<-setNames(data.frame(matrix(ncol = length(columnnames), nrow = 1)), columnnames)

# save the crs of a getGRIDmet product:

GRIDMET_crs<-"GEOGCRS[\"unknown\",\n    DATUM[\"unknown\",\n        ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]],\n    PRIMEM[\"unknown\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433,\n            ID[\"EPSG\",9122]]],\n    CS[ellipsoidal,2],\n        AXIS[\"longitude\",east,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]],\n        AXIS[\"latitude\",north,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]]"

# reproject vect.NWIS:

vect.NWIS.proj<-terra::project(vect.NWIS, GRIDMET_crs)

i<-19

# loop through the watersheds:

for (i in 1:dim(vect.NWIS)[1]){
  # save the samples with their dates for the watershed:
  df.C.j<-filter(df.NWIS.TP_CQ, site_no == vect.NWIS$Name[i])%>%
    drop_na(sample_dt)
  # getGRIDMET only goes back as far as 1979-01-01 so if the sample dates are before that omit site
  # first make sure at least one sample date is good:
  if (sum(df.C.j$sample_dt >= as.Date('1979-01-01'))>=1){
    # set start and end dates for the climate downloads. Add 21 because that is the maximum number of days backto look:
    start_date<-min(df.C.j$sample_dt[df.C.j$sample_dt >= as.Date('1979-01-01')+21])-21
    end_date<-max(df.C.j$sample_dt[df.C.j$sample_dt >= as.Date('1979-01-01')+21])
    # download the climate data for the full range of sample dates for this site:
    rast.pr<-getGridMET(vect.NWIS[i,], varname = vars[1], startDate = start_date, endDate = end_date)[[1]]
    rast.tmax<-getGridMET(vect.NWIS[i,], varname = vars[2], startDate = start_date, endDate = end_date)[[1]]
    rast.tmin<-getGridMET(vect.NWIS[i,], varname = vars[3], startDate = start_date, endDate = end_date)[[1]]
    # if the download doesnt work, omit site: (when i=19, the download function returns and object of class 'glue'. no idea)
    if (class(rast.pr)[1]=="SpatRaster" & class(rast.tmax)[1]=="SpatRaster" & class(rast.tmin)[1]=="SpatRaster"){
      # extract:
      df.pr <- terra::extract(rast.pr, vect.NWIS.proj[i,], mean)
      df.tmax <- terra::extract(rast.tmax, vect.NWIS.proj[i,], mean)
      df.tmin <- terra::extract(rast.tmin, vect.NWIS.proj[i,], mean)
      # create a list of these dfs:
      l.climate<-list(df.pr, df.tmax, df.tmin)
      # reformat to convert dates, and merge with sample dates:
      l.climate<-lapply(seq_along(l.climate), \(i) l.climate[[i]]%>%pivot_longer(cols= starts_with(paste0(vars[i],"_")), names_to = 'Date',values_to = "Value")%>%
                          mutate(Date=as.Date(str_replace(Date, paste0(vars[i],"_"), "")), Climate = vars[i])%>%
                          arrange(desc(Date))%>%
                          select(-ID)
      )
      # create dataframe of these listsand pivot_wider:
      df.climate<-bind_rows(l.climate)%>%
        pivot_wider(names_from = Climate, values_from = Value)
      # determine the sample dates:
      sample_dates<-filter(df.climate, Date %in% df.C.j$sample_dt)
      # create list of the 21 day prior date ranges for the sample dates, and merge the climate variables for the dates for each sample dates 21 day dataframe:
      l.sample_dates<-lapply(sample_dates$Date, \(i) data.frame(Date=seq(i-21, i, by = 1))%>%
                             arrange(desc(Date))%>%
                             left_join(., df.climate, by = 'Date'))
      # perform operations on each sample date using lapply and a function (sourced):
      l.sample_dates<-lapply(l.sample_dates, \(i) fun.Process_climate(i, site = vect.NWIS$Name[i], date = max(i$Date)))
      # bind list of dataframes into a single df:
      df.sample_dates<-bind_rows(l.sample_dates)
      # append this dataframe to master df:
      df.ML_NWIS.Climate_features<-bind_rows(df.ML_NWIS.Climate_features, df.sample_dates)
    }else{print(i); print(vect.NWIS$Name[i]); print(c(class(rast.pr), class(rast.tmax), class(rast.tmin)))}
  }else{print(i); print(vect.NWIS$Name[i]); print('none of site dates will work with getGRIDMet')}
}  
  
# remove the first row
  
df.ML_NWIS.Climate_features<-df.ML_NWIS.Climate_features[-1,]  

# save:

save(df.ML_NWIS.Climate_features, file = 'Processed_Data/df.ML_NWIS.Climate_features.Rdata')

load('Processed_Data/df.ML_NWIS.Climate_features.Rdata')

# takelook around thisdataframe, see which watersheds/dates are not represented:

filter(df.NWIS.TP_CQ, site_no %in% vect.NWIS$Name)%>%distinct(site_no)
unique(df.ML_NWIS.Climate_features$site_no)
# 42 of the 42 sites arerpresented

a<-filter(df.NWIS.TP_CQ, site_no %in% vect.NWIS$Name)%>%group_by(site_no)%>%summarize(n=n())
b<-df.ML_NWIS.Climate_features%>%group_by(site_no)%>%summarize(n=n())

c<-left_join(a,b,by = 'site_no')

sum(c$n.x==c$n.y, na.rm = T) # only two sites worked 'fully' but this is proably due to the samples being before 1979

#### create list of the 5 predictor sets: ####

# read in an edited version of the gauges 2 features (removed time specific climatevariables because those were comming up as highly correlated with TP):

df.G2_edited <- read_excel_allsheets("Raw_Data/gagesII_sept30_2011_conterm_EDITED.xlsx")
df.G2_edited[[24]]<-NULL
df.G2_edited<-reduce(df.G2_edited, full_join, by = "STAID")

# Full predictor set (Non-time + time specific features (sites are in gauges 2 and have delineations, i.e. full amount of data)

Full<-left_join(df.ML_NWIS.Climate_features, df.ML_NWIS.Q_features, by = c('site_no', 'Date'))%>%
  left_join(., df.G2_edited, by = c('site_no'="STAID"))

# non-time specific features + flow (sites are in gauges 2, have flow, but have no delineation)

No_delineations<-left_join(df.ML_NWIS.Q_features, df.G2_edited,  by = c('site_no'="STAID"))

# Just time specific features (sites are not in gauges 2 but have flow and delinations )

Just_time<-left_join(df.ML_NWIS.Climate_features, df.ML_NWIS.Q_features, by = c('site_no', 'Date'))

# Just the flow component of the time specific features (sites are not in gauges 2 and have no delineation, but have flow)

Just_flow<-df.ML_NWIS.Q_features

# Just non-time specific features (sites are in gauges 2 but have no flow and no delineation - I understand that if a site is in gauges 2 it will have flow, but this situation is mimics when we might recreate the gauges 2 features ourselves)

Just_non_time<-filter(df.G2_edited, STAID %in% x$site_no)

# create a list of these dfs:

l.feature_sets<-list(Full, No_delineations, Just_time, Just_flow, Just_non_time)

# rename the date column of raw tp df:

temp<-rename(df.NWIS.TP_CQ, Date = sample_dt)%>%
  select(c(site_no, Date, result_va))

# merge the raw TP data to these dfs: (wouldnt work with lapply, proably causethe conditional needed)

l<-list()

for(i in seq_along(l.feature_sets)){
  print(i)
  if('Date' %in% names(l.feature_sets[[i]])){
    l[[i]]<-left_join(temp, l.feature_sets[[i]], by = c('site_no', 'Date'))
  }
  else {
    l[[i]]<-left_join(temp, l.feature_sets[[i]], by = c('site_no'='STAID'))
  }
  
}


# save.image(file = 'Processed_Data/ML_NWIS.Rdata')

load("Processed_Data/ML_NWIS.Rdata")



#### Correlations for each set #####

# remove non numeric/non date columns:

l.cor<-lapply(l, \(df) df[,-which(sapply(df, class) == "character")])

# remove rows with NAs from each df in list (get complete observations)

l.cor<-lapply(l.cor, \(final) final[complete.cases(final), ])

# filter to just numeric columns:

l.cor<-lapply(l.cor, \(i) i%>%select_if(is.numeric))

# perform cor of dataframe:

l.cor.temp<-lapply(l.cor, \(temp) as.data.frame(cor(temp[-1], temp[1],use="complete.obs", method = 'spearman')))

# create column of row names (which are the predictor names) and set the row names to 1:n:
# arrange by abs(spearman rank correlation):

l.cor.temp.1<-lapply(l.cor.temp, \(temp) temp%>%mutate(Predictor=rownames(temp))%>%arrange(temp, desc(abs(result_va))))


# create a seperate df for the top 35 and add a rankcolumn:
temp<-lapply(l.cor.temp.1, \(temp) temp[1:35,]%>%
  arrange(desc(result_va))%>%
  mutate(Spear_Rank = 1:35)%>%
  mutate(Spear_Rank_name = paste(Spear_Rank, Predictor))
)

# univariate plots:

d1<-temp[[1]]

d1%>%
  select(c(site_no, result_va, temp$Predictor))%>% # select the top 35 predictors (in temp)
  pivot_longer(cols = -c(site_no, result_va), names_to = 'Predictor', values_to = 'Value')%>% # pivot longer to aid ggplot
  mutate(across(Predictor, factor, levels=temp$Predictor))%>% # add the factor levels based on a rank esabilished in temp to get facets to plot in descending spearman rank value order.
  ggplot(., aes(x = Value, y = result_va))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap('Predictor', scales = 'free')

















9# filter the raw Tp CQ df to the 43 sitesIDedabove:

d1<-df.NWIS.TP_CQ%>%filter(site_no %in% x$site_no)

# plot the time series of these sites:

# ggplot(d1, aes(x = sample_dt, y = result_va))+
#   geom_point()+
#   facet_wrap('n_sample_rank')+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank()
#   )

# I would remove the following sites since their time series do nothave
# big chunks of data that WRTDScould be runon:

x<-distinct(d1, site_no, .keep_all = T)

x<-x[-c(25,29,31:34,37,39:43),]

d1<-df.NWIS.TP_CQ%>%filter(site_no %in% x$site_no)

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

# merge the watershed boundaries to this df: to do this:

TP<-left_join(df.sf.NWIS.keep, TP, by = c('Name'='site_no'))%>%
  arrange(desc(cv))%>%
  mutate(cv_rank = 1:nrow(.))

# map:

mapview(TP, zcol = 'cv')

# I want to look at the time series plots ordered by cv:

# add the cv column to the raw tp data df:

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

# perform spearman rank correlations: to do this:

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




save.image(file = 'Processed_Data/ML_NWIS.Rdata')
  
load("Processed_Data/ML_NWIS.Rdata")






