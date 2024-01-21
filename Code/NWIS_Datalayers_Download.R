# Ryan Ruggiero

rm(list=ls(all=T)) # clear global env.
gc()

####################### Load packages #######################

library(FactoMineR)
library(factoextra)
library(leaps)
library(glmnet)
library(randomForest)
library(ggcorrplot)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
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

# download and process watershed attributes from datalayers
# watershed shapefiles come from Code/NWIS_Delineate_and_WWTP.R

####################### Workflow #######################

# read in df.sf.NWIS:

load('Processed_Data/NWIS_Watershed_Shapefiles.Rdata')

####~~~~ Data Layers ~~~~####

#### CDL: ####

# note when I ran ths CDL code block I used the 103 sites thatwere apartofthe set, this is prior to filtering based onthe CDLdate... so I am just going to filter the resulting CDL df at the end. Just note that what is saved is of the 103 sites, but that the code doesnt reflect this, i.e. I just ran the filter anddidnt keep the code because df.NWIS.TP.keep would be down to 56 if I did it right 

# download:

# rast.NWIS.CDL.2020 <- GetCDLData(aoi = df.sf.NWIS, year = "2022", type = "b", tol_time = 1000)
# rast.NWIS.CDL.2008 <- GetCDLData(aoi = df.sf.NWIS, year = "2008", type = "b", tol_time = 1000)

# convert to rast:

# rast.NWIS.CDL.2020<-rast(rast.NWIS.CDL.2020)
# rast.NWIS.CDL.2008<-rast(rast.NWIS.CDL.2008)

# check CRS of both years:

# crs(rast.NWIS.CDL.2020)
# crs(rast.NWIS.CDL.2008)

# plot

# plot(rast.NWIS.CDL.2020)
# plot(rast.NWIS.CDL.2008)

# reproject to sample watershed vector data to match raster data:

# vect.NWIS<-vect(df.sf.NWIS) # need to first ininalize vect since sometimes reading in rdata file (terra issue)
# 
# vect.NWIS.proj<-terra::project(vect.NWIS, crs(rast.NWIS.CDL.2020))

# extract frequency tables for each sample watershed

# l.NWIS.CDL <- terra::extract(rast.NWIS.CDL.2020, vect.NWIS.proj, table,ID=FALSE)[[1]]
# l.NWIS.CDL.2008 <- terra::extract(rast.NWIS.CDL.2008, vect.NWIS.proj, table,ID=FALSE)[[1]]

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
# looking at the CDL legend (CropScapeR::linkdata), I determined the following:

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

# now I am ready to merge withthe legend:
# doing it two ways. The first is the straigth merge where the CDL classes are in acolumn
# and the second way is to use thisreclassified legend: to do this:
# left join each df in the list to the CDL legend key, as well as calcuate the pland:

# the real CDL variables:

l.NWIS.CDL.real<-lapply(l.NWIS.CDL, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., linkdata, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes
l.NWIS.CDL.2008.real<-lapply(l.NWIS.CDL.2008, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., linkdata, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes

# the reclassfied CDL variables (into landuse classes):

l.NWIS.CDL.reclass<-lapply(l.NWIS.CDL, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., reclass_CDL, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes
l.NWIS.CDL.2008.reclass<-lapply(l.NWIS.CDL.2008, \(i) i%>%mutate(Var1 = as.integer(as.character(Var1)),Freq=round(Freq/sum(Freq),2))%>%dplyr::left_join(., reclass_CDL, by = c('Var1' = 'MasterCat'))) # I replaced linkdata in the left join with reclass_CDL to get simplified CDL classes

# set names of list:

names(l.NWIS.CDL.real)<-df.sf.NWIS$Name # note doesnt work with the workflow set up to filter on 2008 limiter
names(l.NWIS.CDL.2008.real)<-df.sf.NWIS$Name # note doesnt work with the workflow set up to filter on 2008 limiter
names(l.NWIS.CDL.reclass)<-df.sf.NWIS$Name # note doesnt work with the workflow set up to filter on 2008 limiter
names(l.NWIS.CDL.2008.reclass)<-df.sf.NWIS$Name # note doesnt work with the workflow set up to filter on 2008 limiter

# combine list into single df:

df.NWIS.CDL<-bind_rows(l.NWIS.CDL.real, .id = 'Name')
df.NWIS.CDL.2008<-bind_rows(l.NWIS.CDL.2008.real, .id = 'Name')
df.NWIS.CDL.reclass<-bind_rows(l.NWIS.CDL.reclass, .id = 'Name')
df.NWIS.CDL.2008.reclass<-bind_rows(l.NWIS.CDL.2008.reclass, .id = 'Name')

# remove potential for one of the MasterCats in the CDL to be empty, which is messing with the pivot_wider below:

df.NWIS.CDL<-filter(df.NWIS.CDL, Crop != '')
df.NWIS.CDL.2008<-filter(df.NWIS.CDL.2008, Crop != '')
df.NWIS.CD.reclassL<-filter(df.NWIS.CDL.reclass, Crop != '')
df.NWIS.CDL.2008.reclass<-filter(df.NWIS.CDL.2008.reclass, Crop != '')

# pivot wider:

# for the CDL linkdata:

df.NWIS.CDL<- pivot_wider(df.NWIS.CDL[,-2], names_from = Crop, values_from = Freq)
df.NWIS.CDL.2008<- pivot_wider(df.NWIS.CDL.2008[,-2], names_from = Crop, values_from = Freq)

# for the reclass_CDL data:

df.NWIS.CDL.reclass<-df.NWIS.CDL.reclass[,-2]%>%
  group_by(Name, Crop)%>%
  summarise(Freq=sum(Freq, na.rm = T))%>%
  pivot_wider(., names_from = Crop, values_from = Freq)

df.NWIS.CDL.2008.reclass<-df.NWIS.CDL.2008.reclass[,-2]%>%
  group_by(Name, Crop)%>%
  summarise(Freq=sum(Freq, na.rm = T))%>%
  pivot_wider(., names_from = Crop, values_from = Freq)

# check to see if add up to 100%:

sort(rowSums(df.NWIS.CDL[,-1], na.rm = T))
sort(rowSums(df.NWIS.CDL.2008[,-1], na.rm = T))
sort(rowSums(df.NWIS.CDL.reclass[,-1], na.rm = T))
sort(rowSums(df.NWIS.CDL.2008.reclass[,-1], na.rm = T))

# potential to find the site if one of the CDL 2020 is low:

which.min(rowSums(df.NWIS.CDL[,-1], na.rm = T))

# comparing 2022 CDL to 2008: todo thisL

# rename dfs:

df1<-df.NWIS.CDL
df2<-df.NWIS.CDL.2008

# replace NAs with zeros:

df1[is.na(df1)]<-0
df2[is.na(df2)]<-0

# remove all zero columns while keeping Name columns:

df1<-df1%>%
  select_if(~!(all(. == 0) && is.numeric(.)))

df2<-df2%>%
  select_if(~!(all(. == 0) && is.numeric(.)))

# use intersect to get the common columns:

cols <- sort(intersect(names(df1), names(df2)))

# remove Name:

cols = cols[!(cols %in% c('Name'))]

# then sort them so that the columns are in same order while subtracting irrespective of their order in their respective data frames

df_temp<-df1[cols] - df2[cols]

# filter down: replace changes that are less than +/- 0.02 to NA first (i.e. those that had 1% change or less I an considering to be nochange)
# also add NAme column back on:

df_temp<-df_temp%>%
  mutate_all(~ifelse(abs(.) < 0.02, NA, .))%>%
  select_if(~any(!is.na(.)))%>%
  mutate(Name = df.NWIS.CDL$Name, .before = 1)
  
# create heatmap of this df: to do this:

# pivot longer:

df_gg<-df_temp%>%pivot_longer(cols = 2:last_col(), values_to = 'Value', names_to = 'Type')

ggplot(df_gg, aes(x = Name, y = Type, fill = Value)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_gradient2(low = "red", high = "green")+
  ggtitle('Heat Map of the CDL Class Changes For 42 NWIS sites between 2008 and 2022\nNegative values represent a loss from 2008 to 2022')

# 01349840 is an outlier site, so removing that one and rerunning plot:

df_gg%>%filter(Name != '01349840')%>%
  ggplot(., aes(x = Name, y = Type, fill = Value)) +
    geom_tile()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    scale_fill_gradient2(low = "red", high = "yellow")+
    ggtitle('Heat Map of the CDL Class Changes For 42 NWIS sites between 2008 and 2022\nNegative values represent a loss from 2008 to 2022')

# that didnt really work, so now ry removing Deciduous_Forest and Mixed_Forest:

df_gg%>%filter(!Type %in% c('Deciduous_Forest','Mixed_Forest'))%>%
  ggplot(., aes(x = Name, y = Type, fill = Value)) +
    geom_tile()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    scale_fill_gradient2(low = "red", high = "green")+
    ggtitle('Heat Map of the CDL Class Changes For 42 NWIS sites between 2008 and 2022\nNegative values represent a loss from 2008 to 2022')


# I think I am done with CDL for now. 
# I am thinking I would use the 2022 CDL in the correlations and MLR 
# but only if the predictors that hadthe most change between 2008 and 2022 
# are in the top correlates/models? Or does it not matter?


# 


















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


# read in watershed shapefiles:

load('Processed_Data/NWIS_Watershed_Shapefiles.Rdata')







