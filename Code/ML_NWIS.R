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

filter(df.sf.NWIS.keep, Name %in% x$site_no)%>%mapview(.)

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

# loop through the climate variables:

i<-19
j<-1
k<-1

# loop through the watersheds:

for (i in 1:dim(vect.NWIS)[1]){
  # save the samples with their dates for the watershed:
  df.C.j<-filter(df.NWIS.TP_CQ, site_no == vect.NWIS$Name[i])%>%
    drop_na(sample_dt)
  # download the climate data for the full range of sample dates for this site:
  rast.pr<-getGridMET(vect.NWIS[i,], varname = vars[1], startDate = as.character(min(df.C.j$sample_dt[j])), endDate = as.character(max(df.C.j$sample_dt[j])))[[1]]
  rast.tmax<-getGridMET(vect.NWIS[i,], varname = vars[1], startDate = as.character(min(df.C.j$sample_dt[j])), endDate = as.character(max(df.C.j$sample_dt[j])))[[1]]
  rast.tmin<-getGridMET(vect.NWIS[i,], varname = vars[1], startDate = as.character(min(df.C.j$sample_dt[j])), endDate = as.character(max(df.C.j$sample_dt[j])))[[1]]
  
  # loop through the sample dates for the watershed:
  for (j in seq_along(df.C.j$sample_dt)){
    # loop through the climate variables
    for (k in seq_along(vars)){
      # getGRIDMET only goes back as far as 1979-01-01 so if the sample date is before that omit:
      if (df.C.j$sample_dt[j] >= as.Date('1979-01-01')){
        # download:
        rast.climate.k<-getGridMET(vect.NWIS[i,], varname = vars[k], startDate = as.character(df.C.j$sample_dt[j]-21), endDate = as.character(df.C.j$sample_dt[j]))[[1]]
        # if the download doesnt work, omit: (when i=19, the download function returns and object of class 'glue'. no idea)
        if (class(rast.climate.k)[1]=="SpatRaster"){
          # reproject for extract:
          vect.NWIS.proj<-terra::project(vect.NWIS[i,], crs(rast.climate.k))
          # extract:
          df.climate.k <- terra::extract( rast.climate.k, vect.NWIS.proj, mean)
          # reformat to convert dates:
          df.climate.k.1<-df.climate.k%>%
            pivot_longer(cols= starts_with(paste0(vars[k],"_")), names_to = 'Date',values_to = "Value")%>%
            mutate(Date=as.Date(str_replace(Date, paste0(vars[k],"_"), "")))%>%
            arrange(desc(Date))
          # use logical statements to determine what to do this with info:
          if (k == 1){
            # create dataframe of the rain column names and theirvalues using cumsum and subseting vector using the vector of rain_cum_and_temp_deltas (prior to loo)
            df.rain<-data.frame(Name = rain_colnames, Value = round(cumsum(df.climate.k.1$Value)[rain_cum_and_temp_deltas],2))%>%
              pivot_wider(names_from = Name, values_from = Value)%>%
              mutate(site_no = vect.NWIS$Name[i], Date = df.C.j$sample_dt[j], .before = 1)
          }
          if (k == 2){
            df.tmax.lag<-data.frame(Name = tmax_lag_colnames, Value = df.climate.k.1$Value[temp_lags+1])%>%
              pivot_wider(names_from = Name, values_from = Value)%>%
              mutate(site_no = vect.NWIS$Name[i], Date = df.C.j$sample_dt[j], .before = 1)
            df.tmax.delta<-data.frame(Name = tmax_delta_colnames, Value = round(df.climate.k.1$Value[rain_cum_and_temp_deltas+1]-df.climate.k.1$Value[1],2))%>%
              pivot_wider(names_from = Name, values_from = Value)%>%
              mutate(site_no = vect.NWIS$Name[i], Date = df.C.j$sample_dt[j], .before = 1)
            df.tmax<-merge(df.tmax.lag, df.tmax.delta)
          }
          if (k == 3){
            df.tmin.lag<-data.frame(Name = tmin_lag_colnames, Value = df.climate.k.1$Value[temp_lags+1])%>%
              pivot_wider(names_from = Name, values_from = Value)%>%
              mutate(site_no = vect.NWIS$Name[i], Date = df.C.j$sample_dt[j], .before = 1)
            df.tmin.delta<-data.frame(Name = tmin_delta_colnames, Value = df.climate.k.1$Value[rain_cum_and_temp_deltas+1]-df.climate.k.1$Value[1])%>%
              pivot_wider(names_from = Name, values_from = Value)%>%
              mutate(site_no = vect.NWIS$Name[i], Date = df.C.j$sample_dt[j], .before = 1)
            df.tmin<-merge(df.tmin.lag, df.tmin.delta)
          }
        } 
        else {
          # print(vect.NWIS$Name[i])
          # print(df.C.j$sample_dt[j])
          # print(vars[k])
          # print('bad download')
          df.rain<-structure(list(site_no = NA_character_, Date = structure(NA_real_, class = "Date"), 
                                  Rain_cummulative_1_day = NA_real_, Rain_cummulative_2_day = NA_real_, 
                                  Rain_cummulative_3_day = NA_real_, Rain_cummulative_4_day = NA_real_, 
                                  Rain_cummulative_5_day = NA_real_, Rain_cummulative_10_day = NA_real_, 
                                  Rain_cummulative_14_day = NA_real_, Rain_cummulative_21_day = NA_real_), row.names = c(NA, 
                                                                                                                         -1L), class = c("tbl_df", "tbl", "data.frame"))
          df.tmax<-structure(list(site_no = NA_character_, Date = structure(NA_real_, class = "Date"), 
                                  Tmax_lag_0_day = NA_real_, Tmax_lag_1_day = NA_real_, Tmax_lag_2_day = NA_real_, 
                                  Tmax_lag_3_day = NA_real_, Tmax_lag_4_day = NA_real_, Tmax_lag_5_day = NA_real_, 
                                  Tmax_lag_10_day = NA_real_, Tmax_lag_14_day = NA_real_, Tmax_lag_21_day = NA_real_, 
                                  Tmax_delta_1_day = NA_real_, Tmax_delta_2_day = NA_real_, 
                                  Tmax_delta_3_day = NA_real_, Tmax_delta_4_day = NA_real_, 
                                  Tmax_delta_5_day = NA_real_, Tmax_delta_10_day = NA_real_, 
                                  Tmax_delta_14_day = NA_real_, Tmax_delta_21_day = NA_real_), row.names = c(NA, 
                                                                                                             -1L), class = "data.frame")
          df.tmin<-structure(list(site_no = NA_character_, Date = structure(NA_real_, class = "Date"), 
                                  Tmin_lag_0_day = NA_real_, Tmin_lag_1_day = NA_real_, Tmin_lag_2_day = NA_real_, 
                                  Tmin_lag_3_day = NA_real_, Tmin_lag_4_day = NA_real_, Tmin_lag_5_day = NA_real_, 
                                  Tmin_lag_10_day = NA_real_, Tmin_lag_14_day = NA_real_, Tmin_lag_21_day = NA_real_, 
                                  Tmin_delta_1_day = NA_real_, Tmin_delta_2_day = NA_real_, 
                                  Tmin_delta_3_day = NA_real_, Tmin_delta_4_day = NA_real_, 
                                  Tmin_delta_5_day = NA_real_, Tmin_delta_10_day = NA_real_, 
                                  Tmin_delta_14_day = NA_real_, Tmin_delta_21_day = NA_real_), row.names = c(NA, 
                                                                                                             -1L), class = "data.frame")
        }
      }
      else{
        # print(vect.NWIS$Name[i])
        # print(df.C.j$sample_dt[j])
        # print(vars[k])
        # print('bad date')
        df.rain<-structure(list(site_no = NA_character_, Date = structure(NA_real_, class = "Date"), 
                                Rain_cummulative_1_day = NA_real_, Rain_cummulative_2_day = NA_real_, 
                                Rain_cummulative_3_day = NA_real_, Rain_cummulative_4_day = NA_real_, 
                                Rain_cummulative_5_day = NA_real_, Rain_cummulative_10_day = NA_real_, 
                                Rain_cummulative_14_day = NA_real_, Rain_cummulative_21_day = NA_real_), row.names = c(NA, 
                                                                                                                       -1L), class = c("tbl_df", "tbl", "data.frame"))
        df.tmax<-structure(list(site_no = NA_character_, Date = structure(NA_real_, class = "Date"), 
                                Tmax_lag_0_day = NA_real_, Tmax_lag_1_day = NA_real_, Tmax_lag_2_day = NA_real_, 
                                Tmax_lag_3_day = NA_real_, Tmax_lag_4_day = NA_real_, Tmax_lag_5_day = NA_real_, 
                                Tmax_lag_10_day = NA_real_, Tmax_lag_14_day = NA_real_, Tmax_lag_21_day = NA_real_, 
                                Tmax_delta_1_day = NA_real_, Tmax_delta_2_day = NA_real_, 
                                Tmax_delta_3_day = NA_real_, Tmax_delta_4_day = NA_real_, 
                                Tmax_delta_5_day = NA_real_, Tmax_delta_10_day = NA_real_, 
                                Tmax_delta_14_day = NA_real_, Tmax_delta_21_day = NA_real_), row.names = c(NA, 
                                                                                                           -1L), class = "data.frame")
        df.tmin<-structure(list(site_no = NA_character_, Date = structure(NA_real_, class = "Date"), 
                                Tmin_lag_0_day = NA_real_, Tmin_lag_1_day = NA_real_, Tmin_lag_2_day = NA_real_, 
                                Tmin_lag_3_day = NA_real_, Tmin_lag_4_day = NA_real_, Tmin_lag_5_day = NA_real_, 
                                Tmin_lag_10_day = NA_real_, Tmin_lag_14_day = NA_real_, Tmin_lag_21_day = NA_real_, 
                                Tmin_delta_1_day = NA_real_, Tmin_delta_2_day = NA_real_, 
                                Tmin_delta_3_day = NA_real_, Tmin_delta_4_day = NA_real_, 
                                Tmin_delta_5_day = NA_real_, Tmin_delta_10_day = NA_real_, 
                                Tmin_delta_14_day = NA_real_, Tmin_delta_21_day = NA_real_), row.names = c(NA, 
                                                                                                           -1L), class = "data.frame")
      }
    }
    # combine the three dataframes made in the k prior loop:
    df.ij.row<-merge(df.rain, df.tmax)%>%merge(.,df.tmin)
    # append to master df:
    df.ML_NWIS.Climate_features<-bind_rows(df.ML_NWIS.Climate_features, df.ij.row)
  }
}

# didnt finish, I stopped it: i @37, j@989

# save what I have so far:

# save(df.ML_NWIS.Climate_features, file = 'Processed_Data/df.ML_NWIS.Climate_features.Rdata')

load('Processed_Data/df.ML_NWIS.Climate_features.Rdata')

# take look around this dataframe, see which watersheds/dates are not represented:



# create list of the 5 predictor sets:

Full<-


















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






