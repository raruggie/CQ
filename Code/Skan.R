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

# Step1: delinate watersheds
# Step 2: download surrogate gauge flow data (Owasco Inlet)
# Step3: import C data to see what data from the data layers are needed 
# Step4: collect watershedcharacterstics from data layers
# Step 5: combine C, Q, and data layers into a single data frame
# Step 6: fit CQ curves to sites
# Step 7: perform correlations with CQ parameters and watershed characteristics

# Side goal: data frame is ready for ML and code is flexible to add sites and sample dates

####################### Workflow #######################

# import site locations:

df.Skan.TP.sites <- read_excel("C:/PhD/CQ/Raw_Data/AISkanListTP.xlsx")

# filter the sites to those with lat longs and combine trib name and type into 23 unique names:

df.Skan.TP.sites<-df.Skan.TP.sites%>%drop_na(Latitude)%>%
  mutate(`Site Name`=paste(`Site Name`,Type))

# delinate each site and run through steps to convert DA geometries. todothis: 

# use function (source) todownload DA and set the names of the DA sites in the list:

l.SS_WS.Skan<-lapply(seq_along(df.Skan.TP.sites$Latitude), \(i) Ryan_delineate(df.Skan.TP.sites$Longitude[i],df.Skan.TP.sites$Latitude[i]))

names(l.SS_WS.Skan)<-df.Skan.TP.sites$`Site Name`

save(l.SS_WS.Skan, file = 'C:/PhD/CQ/Downloaded_Data/l.SS_WS.Skan.Rdata')

# convert to df sf using function (source):

df.sf.Skan<-fun.l.SS_WS.to.sfdf(l.SS_WS.Skan)

save(df.sf.Skan, file = 'C:/PhD/CQ/Processed_Data/df.sf.Skan.Rdata')

# look at a map of the watersheds:

mapview(df.sf.Skan, zcol = 'Name')

# add column with drainage area:
# now you are working in vect 

vect.Skan<-vect(df.sf.Skan)

vect.Skan$area_KM2<-expanse(vect.Skan, unit="km")

# look at vect df:

head(vect.Skan[,c('Name', 'area_KM2')], 25)

# want to look at the stream network:

# read in Skan lake watershed boundary:

df.sf.Skan_DA<-st_read('C:/PhD/Research/Data/Skan/Skan_watershed_shapefile/globalwatershed.shp')

# look at map of basin:

mapview(df.sf.Skan_DA)

# now use this boundary to get national hydrograph dataset (includes stream network): to dothis:

# first need to make valid:

df.sf.Skan_DA<-st_make_valid(df.sf.Skan_DA)

#download nhd:

# nhd.Skan<-nhdget_nhd(df.sf.Skan_DA,label = '2') # just takes a while to run

# save(nhd.Skan, file='C:/PhD/CQ/Downloaded_Data/nhd.Skan.Rdata')

# save justthe flow path:

nhd.Skan.fp<-nhd.Skan$Flowline

# filteroutsmallerstreams? notsure ifthis works so left commented out (not in saved file):
# nhd_fp<-filter(nhd_fp, !is.na(gnis_name))

# look at nhd map

mapview(df.sf.Skan,zcol = 'Name')+mapview(nhd.Skan.fp)

####################### Part 2 - Import flow data and scale to each site #######################

# using Owasco inlet to scale flows:

# download the draiange area size from site metadata
# then convert draiange area to km2

df.Ow<-readNWISsite('04235299')

v.Ow_DA<-df.Ow$drain_area_va*2.58999

# download raw flow data:

df.Ow_Q<-readNWISdv('04235299', '00060')

# convert date to Date and filter the data from 1960s:

df.Ow_Q<-df.Ow_Q%>%mutate(Date = as.Date(Date))%>%
  filter(Date > as.Date('2000-01-01'))%>%
  select(c(3,4))%>%
  rename(Q = 2)

# I will transfer this flow into a list of multiple scaled flow dfs later

####################### Part 3 - Import C data to see what is needed for data layers #######################

# I will do this later when I get the raw TP data

#### Step 1: Import C data

df.Skan.C <- read_excel("C:/PhD/Research/Data/Skan/SkanStream-Ryan_sites_edited_to_match_AISkanListTP.xlsx")

# what are the min and max dates for CDL, climate data layers?

min(df.Skan.C$Date) # 2018

max(df.Skan.C$Date) # 2021

# this is a small range for CDL, so I am just going to use 1 CDL year instead of 4
# but climate data will be every year

####################### Part 4 - Data Layers Workflow #######################

# What we want:
# a dataframe with rows for each site and columns for each watershed characteristic

# We start with downloading SpatRasters of the entire Skan draiange area
# I was orginally thinking merging them all in a single SPatRaster object as layers,
# but they have different CRS, so that doesnt make sense.

# So, one option is to create a workflow that crops the raster one at a time for each subbwatershed (as entires in a SpatVectors)
# then run the calculations on the cropped rasters for each watershed 

# so the workflow for each raster is:
# 1) download and convert to SpatRast if not already
# 2) extract data from the raster at each polygon location

#### CDL :

# download:

rast.Skan.CDL.2020 <- GetCDLData(aoi = df.sf.Skan_DA, year = "2020", type = "b", tol_time = 1000)

# convert to rast:

rast.Skan.CDL.2020<-rast(rast.Skan.CDL.2020)

# plot

plot(rast.Skan.CDL.2020)

# extract frequency tables for each sample watershed
# first reproject to sample watershed vector data to match raster data:
# then extract raster frequencies over polygons
# then convert resulting list of tables to list of dfs
# and finally left join each df in the list to the CDL legend key, as well as calcuate the pland:

vect.Skan.proj<-terra::project(vect.Skan, crs(rast.Skan.CDL.2020))

l.Skan.CDL <- terra::extract(rast.Skan.CDL.2020, vect.Skan.proj, table,ID=FALSE)[[1]]

l.Skan.CDL<-lapply(l.Skan.CDL, as.data.frame)

l.Skan.CDL<-lapply(l.Skan.CDL, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., linkdata, by = c('Var1' = 'MasterCat')))

names(l.Skan.CDL)<-df.sf.Skan$Name

df.Skan.CDL<-bind_rows(l.Skan.CDL, .id = 'Name')

df.Skan.CDL<-df.Skan.CDL%>%
  select(!Var1)%>%
  pivot_wider(., names_from = Crop, values_from = Freq)

#### NED

DEM<-get_ned(df.sf.Skan_DA, label = '2') # already SpatRaster!

plot(DEM)

# extract elevation metrics over each sample watershed:
# can use a function with multiple functions. If not doing this, you would need to run extract for each metric:

vect.Skan.proj<-terra::project(vect.Skan, crs(DEM))

f <- function(x, na.rm = T) {
  c(mean=mean(x, na.rm = na.rm),
    range=max(x, na.rm = na.rm)-min(x, na.rm = na.rm),
    sd=sd(x, na.rm = na.rm)
  )
}

df.Skan.DEM <- as.data.frame(terra::extract(DEM, vect.Skan.proj, f))

names(df.Skan.DEM)<-c('Name', 'Elev_Avg', 'Elev_Range', 'Elev_SD')

df.Skan.DEM$Name<-df.sf.Skan$Name

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

df.Skan.Climate<-data.frame(ID = 1:dim(df.Skan.CDL)[1])

i<-1

for (i in seq_along(vars)){
  
  # download:
  
  climate<-getGridMET(vect.Skan, varname = vars[i], startDate = "2016-01-01", endDate = "2019-12-31")[[1]]
  
  # extract:
  
  vect.Skan.proj<-terra::project(vect.Skan, crs(climate))
  
  climate <- terra::extract(climate, vect.Skan.proj, mean)
  
  # reformat, convert dates, and calcualte annual stats:
  
  climate1<-climate%>%
    pivot_longer(cols= starts_with(paste0(vars[i],"_")), names_to = 'Date',values_to = "Value")%>%
    mutate(Date=as.Date(str_replace(Date, paste0(vars[i],"_"), "")))%>%
    group_by(year(Date),ID)%>%
    summarize(Annual = vars_funs[[i]](Value))%>%
    rename(Year=1)%>%
    mutate(Year=paste0(vars[i],Year))%>%
    pivot_wider(names_from = Year, values_from = Annual)
  
  df.Skan.Climate<-left_join(df.Skan.Climate, climate1, by = 'ID')
  
}

# rename the ID column

df.Skan.Climate[,1]<-df.Skan.DEM$Name

names(df.Skan.Climate)[1]<-'Name'

##### Finally combine CDL, DEM, climate, and Streamstats predictors into common df:

df.Skan.Predictors<-left_join(df.Skan.CDL,df.Skan.DEM, by = 'Name')%>%left_join(.,df.Skan.Climate, by = 'Name')%>%left_join(.,df.sf.Skan, by = 'Name')

save(df.Skan.Predictors, file = "C:/PhD/CQ/Processed_Data/df.Skan.Predictors.Rdata")

####################### Part 5 - combine C,Q, and data layers #######################

# Need to combine C data set site format with data layers dataframe 
# first reduce df.Skan.Predictors to just unique sites
# I cofirmed that the sites that have both SU and UFI names have similar datalayer entries by visual examinaiton of the dataframe

# create a new column with the SU/UFI removed, then create a distinct dataframe from the name column:

df.Skan.Predictors<-df.Skan.Predictors%>%
  mutate(Name= trimws(substr(Name,1,nchar(Name)-3)), .after = 1)%>%
  distinct(Name, .keep_all = T)

# Note that I cleaned up the df.Skan.C data in excel and saved as a new file, reading it in again here:

df.Skan.C <- read_excel("C:/PhD/Research/Data/Skan/SkanStream-Ryan_sites_edited_to_match_AISkanListTP.xlsx")

# lets check the flows: there is a column with the Owasco INlet flows in df.Skan.C
# but I want to check to see if they have been DA scaled:
# convert cfs to cms:

df.Skan.C<-df.Skan.C%>%
  mutate(Date = as.Date(Date))%>%
  left_join(.,df.Ow_Q, by = 'Date')%>%
  relocate(Q, .after=4)%>%
  rename(Name = Location)%>%
  mutate(Q_cms = Q*0.028316832, .after=5 )

# they are the same so not DA scaled
# to DA scale, need to add the drainage areas and drainage area ratios (Skan subbasin DA/Ow DA) to data layers:

vect.Skan<-vect(df.sf.Skan)

vect.Skan$area_KM2<-expanse(vect.Skan, unit="km")

df.Skan.Predictors<-as.data.frame(vect.Skan[,c('Name', 'area_KM2')])%>%
  mutate(Name= trimws(substr(Name,1,nchar(Name)-3)), .after = 1)%>%
  distinct(Name, .keep_all = T)%>%
  left_join(.,df.Skan.Predictors, by = 'Name')%>%
  mutate(sub_Ow_DA_ratio = area_KM2/v.Ow_DA, .after = 2)

# now merge the area ratios to df.Skan.C andcalcualte drainage area scaled flow:

df.Skan.C<-left_join(df.Skan.Predictors[,c('Name','area_KM2',"sub_Ow_DA_ratio")], df.Skan.C, by = 'Name')%>%
  mutate(Q_cms_scaled = Q_cms*sub_Ow_DA_ratio, .after = 4)%>%
  select(!c(2,3,8,9)) # remove duplicated columns in data layers df and unwanted flow columns

# now merge the data layers with the CQ dataframe

df.Skan.CQ.with.Predictors<-left_join(df.Skan.C, df.Skan.Predictors, by = 'Name')

# note that if I want antcdent conditions I will need to very much changethe climate for loop (proably a double for loop)
# note that the number of unique sites here is not as many as in the first few parts of the code because not all the sampling sites lined up with the sites in AISKanTPlist.csv

save.image(file = "C:/PhD/CQ/Processed_Data/Skan.Rdata")

# saved up to here
























####################### Part 6 - Fit 2 slope models to TP CQ curves #######################

# look at CQ curves for TP:

df.Part6<-df.CQ_DL%>%filter(!is.na(`Total Phosphorus (TP) umol/L (From ICP-MS)`))

ggplot(data = df.Part6, aes(x = log(Q_cms_scaled), y = log(`Total Phosphorus (TP) umol/L (From ICP-MS)`)))+
  geom_smooth(method = 'lm')+
  geom_point()+
  facet_wrap('Name', scales = 'fixed')

# fit CQ curves to each site:

# does a two slope model better fit the data?
# this is the alternative hypothesis of the davies test conducted below

# create linear model for each site using the list of sites made above
# then run breakpoint analysis (segmented function) and davies test on the lm's
# create a dataframe for the loop iteration and append to new list
# also, the davies test results are saved to a matrix
# first need to modify the lists and use a for loop - segmeneted function is not liking lapply approach

l_Seg<-list()

davies.test.matrix<-NULL

# i<-unique(df.Part6$Name)[4]

for (i in sort(unique(df.Part6$Name))){
  print(i)
  
  # create a dataframe that will work with segmented
  df<-df.Part6%>%
    filter(Name == i)%>%
    mutate(log_C = log(`Total Phosphorus (TP) umol/L (From ICP-MS)`), log_Q = log(Q_cms_scaled), C = `Total Phosphorus (TP) umol/L (From ICP-MS)`, Q = Q_cms_scaled)%>%
    filter(is.finite(log_C))
  
  # build lm
  m<-lm(log_C~log_Q, df)
  
  # perform breakpoint regression
  m_seg<-segmented(obj = m, npsi = 1)
  
  # and test for constant linear predictor:
  
  x<-paste(i, '-', davies.test(m)$p.val<0.05)
  
  davies.test.matrix<-c(davies.test.matrix,x)
  
  # get the breakpoints
  bp<-m_seg$psi[1]
  
  # get the slopes

  if(length(class(m_seg))==2){
    s<-as.data.frame(slope(m_seg))
  } else{
    s<-NA
  }
  
  # get the intercepts
  
  if(length(class(m_seg))==2){
    inter<-as.data.frame(intercept(m_seg))
  } else{
    inter<-NA
  }
  
  
  # get the fitted data
  
  fit <- data.frame(Q = df$log_Q, Seg_C = fitted(m_seg))
  
  # create a dataframe of the data to export out of the loop
  
  if(length(class(m_seg))==2){
    result_df<-fit%>%mutate(site = i, Q_real = df$Q, C = df$C, I1 = inter$Est.[1], I2 = inter$Est.[2], Slope1 = s[1,1], Slope2 = s[2,1], BP = bp)
  } else{
    result_df<-fit%>%mutate(site = i, Q_real = df$Q, C = df$C, I1 = NA, I2 = NA, Slope1 = NA, Slope2 = NA, BP = NA)
  }
  
  l_Seg[[i]]<-result_df
  
}

# look at the davies test result:

davies.test.matrix

# make a plot of the results
# first combine the list df into a single df,
# then combine the C-Q data and the model fit data into a new 'all' dataframe
# then create a dataframe to plot text
# finally make plot:

df_Seg<-bind_rows(l_Seg)

BP1_p_vals<-trimws(gsub(".*-","",davies.test.matrix))

BP1_text <- data.frame(Q_real= 0.082085,C = 0.082085,lab = BP1_p_vals,
                       site = factor(sort(unique(df.Part6$Name)),levels = sort(unique(df.Part6$Name))))


ggplot(df_Seg, aes(x = log(Q_real), y = log(C)))+
  geom_smooth(method = 'lm')+
  geom_point()+
  geom_line(aes(x = Q, y = Seg_C), color = 'tomato')+
  facet_wrap(dplyr::vars(site))+
  geom_text(data = BP1_text,label = BP1_p_vals)

# export the BP analysis:
# add a logical column to the dataframe for if there are two slopes
# need to first transform the davies.test.matrix into a dataframe to use as a key

dav_df<-as.data.frame(davies.test.matrix)

dav_df<-as.data.frame(str_split_fixed(dav_df$davies.test.matrix, " - ", 2))

dav_df<-dav_df%>%rename(site = 1, BP_yes = 2)

df_Seg<-left_join(df_Seg, dav_df, by = 'site')

df_Seg_dis<-distinct(df_Seg, site, .keep_all = T)%>%select(., -c(1,2,4,5))

# for the sites that dont need two slope models, fit single slope models:

smdf<-df_Seg_dis%>%filter(BP_yes == 'FALSE') # single model df

# now run linear models for each site ad extract intercept andslopes
smdf<-df.Part6%>%
  filter(Name %in% smdf$site)%>%
  mutate(log_C = log(`Total Phosphorus (TP) umol/L (From ICP-MS)`), log_Q = log(Q_cms_scaled), C = `Total Phosphorus (TP) umol/L (From ICP-MS)`, Q = Q_cms_scaled)%>%
  filter(is.finite(log_C))%>%
  group_by(Name)%>%
  do({ co <- coef(lm(log_C ~ log_Q, .))
  summarize(., I = co[1], 
            S = co[2])
  }) %>%
  ungroup

# then merge the I and S back to df_Seg_dis:

df_Seg_dis<-left_join(df_Seg_dis, smdf, by = c('site'='Name'))

####################### Part 7 - Correlations #######################

# Given that there is a mix between two and one slope CQ curves,
# correlaitons will be looked at with just single slope models
# Two models will be used, OLS and Sens slope:

# create a dataframe with OLS and sens slope intercept and slopes:

df.TP.CQ_OLS<-df.Part6%>%
  mutate(log_C = log(`Total Phosphorus (TP) umol/L (From ICP-MS)`), log_Q = log(Q_cms_scaled), C = `Total Phosphorus (TP) umol/L (From ICP-MS)`, Q = Q_cms_scaled)%>%
  filter(is.finite(log_C))%>%
  group_by(Name)%>%
  do({ OLS.co <- coef(lm(log_C ~ log_Q, .))
  summarize(., OLS.I = OLS.co[1], 
            OLS.S = OLS.co[2])
  }) %>%
  ungroup

df.TP.CQ_Sens<-df.Part6%>%
  mutate(log_C = log(`Total Phosphorus (TP) umol/L (From ICP-MS)`), log_Q = log(Q_cms_scaled), C = `Total Phosphorus (TP) umol/L (From ICP-MS)`, Q = Q_cms_scaled)%>%
  filter(is.finite(log_C))%>%
  group_by(Name)%>%
  do({ Sens.co<-zyp.sen(log_C~log_Q,.)
  summarize(., Sen.I = Sens.co$coefficients[[1]],
            Sen.S= Sens.co$coefficients[[2]])
  }) %>%
  ungroup

# merge OLS and Sens:

df.TP.CQ_slopes<-left_join(df.TP.CQ_OLS, df.TP.CQ_Sens, by = 'Name')

# now we can merge the watershed characteristic data to this dataframe

df.TP.CQ_slopes<-left_join(df.TP.CQ_slopes, df.Skan.Predictors, by = 'Name')

# now run correlations between intercepts and slopes and watershed characteristics

# create a dataframe with just the preditor and predictands 
# create variable with the number of months

# going to need to do this four times, one for each intercept andslope for OLS and Sens:

df.TP.cor<-df.TP.CQ_slopes #[,-c(1:4)]

n_month<-dim(df.TP.cor)[1] # just carrying over the name n_months from where I got this code, C:\PhD\Research\Mohawk\Code\Mohawk_Regression-analyizing_predictor_df.R

# now use the corrr package to correlate() and focus() on your variable of choice

df.TP.cor <- df.TP.cor %>% 
  correlate() %>% 
  focus(c(OLS.I, OLS.S, Sen.I, Sen.S))%>%
  pivot_longer(cols= c(2:5), names_to = 'CQ_Parameter', values_to = 'Pearson_Correlation')%>%
  mutate(p_val = round(2*pt(-abs(Pearson_Correlation*sqrt((n_month-2)/(1-(Pearson_Correlation)^2))), n_month-2),2))%>%
  mutate(sig_0.05 = ifelse(p_val <= 0.05, 'sig', 'not'))%>%
  drop_na(p_val) # some standard deviaitons return NA because the watershed characteristic values are zero

# then ggplot2 package to plot the results

l.cor<-df.TP.cor %>%
  split(., df.TP.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>%mutate(term = factor(term, levels = unique(term[order(Pearson_Correlation)])))%>%filter(!between(Pearson_Correlation, -0.25,.25))) # Order by correlation strength
  
plist<-lapply(l.cor, \(i) i%>%ggplot(aes(x = term, y = Pearson_Correlation, color = sig_0.05)) +
         geom_bar(stat = "identity") +
         facet_wrap('CQ_Parameter')+
         ylab(paste('Pearson Correlation')) +
         xlab("Watershed Attribute")+
         theme(axis.text.x=element_text(angle=40,hjust=1))+
         theme(legend.position="bottom"))

ggpubr::ggarrange(plotlist = plist, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")





load("C:/PhD/CQ/Processed_Data/Skan_CQ.Rdata")




