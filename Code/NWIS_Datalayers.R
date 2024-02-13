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

# 1) recreate the GAGES II predictors using data layers for the 
# 42 sites I ID in the NWIS workflow, using the GAGES values for QA

# 2) (actually comes first) comparing early and recent CDL values

####################### Workflow #######################

# read in watershed shapsfiles (df.sf.NWIS):

load('Processed_Data/NWIS_Watershed_Shapefiles.Rdata')

#### CDL Comparisons: Early and Recent downloads and comparisons ~~~####

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

# read in CDL for 2022 and 2008 (recent and old):

load('Processed_Data/l.NWIS.CDL.Rdata')
load('Processed_Data/l.NWIS.CDL.2008.Rdata')

# convert resulting list of tables to list of dfs

l.NWIS.CDL<-lapply(l.NWIS.CDL, as.data.frame)
l.NWIS.CDL.2008<-lapply(l.NWIS.CDL.2008, as.data.frame)

# the next step is to join the CDL key with these dfs, but first:

x<-CropScapeR::linkdata

# aggregate CDL: to do this:
# looking at the CDL legend (CropScapeR::linkdata), I determined the following:

Ag<-c(1:6,10:14,21:39,41:61,66:72,74:77,204:214,216:227,229:250,254, 176, 92) # I added Pasture to this to get it more similar to PLANTNLCD06 from GAGES II
Pasture<-c(36,37,58,176) # I added some other pasture sounding ones to get is closer to PASTURENLCD06
Forest<-c(63,141:143)
Developed<-c(82,121:124)
Water<-c(83,111)
Wetlands_all<-c(87,190,195)
Other<-c(64,65,88,112,131,152) # Shrub<-c(64,152) Barren<-c()

l <- tibble::lst(Ag,Pasture,Forest,Developed,Water,Wetlands_all,Other)

reclass_CDL<-data.frame(lapply(l, `length<-`, max(lengths(l))))%>%
  pivot_longer(cols = everything(), values_to = 'MasterCat',names_to = 'Crop')%>%
  drop_na(MasterCat)

# now I am ready to merge with the legend:
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
df.NWIS.CDL.reclass<-filter(df.NWIS.CDL.reclass, Crop != '')
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

# pivot longer:

# create heatmap of this df:

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

# that didnt really work, so now try removing Deciduous_Forest and Mixed_Forest:

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


























#### Make Table of Predictors used in Manuscript ####

# these are the GAES II predictors:

load('Processed_Data/df.G2.reduced.Rdata')

v<-names(df.G2.reduced)

# HYDRO_DISTURB_INDX can be replaced by FRAGUN_BASIN
# pre1990_STOR may not be able to get but going to try using NHD
# HIRES_LENTIC_NUM, HIRES_LENTIC_DENS, and HIRES_LENTIC_MEANSIZ (number,density, and mean size of lakes ponds resvioirs, maybe just need 1?) can get from NHD
# DEVNLCD06 = Developed, FORESTNLCD06 = Forest, PLANTNLCD06 = Ag+Pasture, WATERNLCD06 = Water, SNOWICENLCD06+BARRENNLCD06+GRASSNLCD06 = Other, DECIDNLCD06+EVERGRNLCD06+MIXEDFORNLCD06+SHRUBNLCD06 = Forest, PASTURENLCD06 = Pasture,CROPSNLCD06 = Ag,WOODYWETNLCD06+EMERGWETNLCD06=Wetlands_All
# CDL_CORN,CDL_SOYBEANS,CDL_SPRING_WHEAT,CDL_WINTER_WHEAT,CDL_WWHT_SOY_DBL_CROP,CDL_ALFALFA,CDL_OTHER_HAYS, are covered in CDL_ALL_OTHER_LAND 

# make a table of predictors: to do this:

# read in variable description:

var_desc <- readxl::read_excel("Raw_Data/gagesII_sept30_2011_var_desc.xlsx")

# filter to those in df.G2.reduced create counter column, and reduce columns:

var_desc<-filter(var_desc, VARIABLE_NAME %in% v)%>%
  mutate(Count = row_number())%>%
  select(Count, VARIABLE_NAME, DESCRIPTION)

# make kable table:

var_desc %>%
  kbl(align = "c",escape = F) %>%
  kable_classic(html_font = 'Times', font_size = 14, full_width = F)%>%
  row_spec(seq(1,nrow(var_desc),2), background="#FF000020")

#





























#### Test of howclose I can get CDL to NLCD ####

# **NOTE*
# the CDL and NLCD differ in that the CDL combined Pasture and grassland while the nLCD does not
# thus, even though the CDL can be agrgated to get close to the NLCD, I would need to do NLCD to get exactly like GAGES
# here is a test of howclose I can get CDL to NLCD

# the CDL reclass are in this df:

df.NWIS.CDL.2008.reclass

# compare them to the GAGES II NLCD ones (again, where we can, see above note):to do this:

temp<-df.G2.reduced%>%select(STAID, names(df.G2.reduced)[c(8:11,18,19,21,22)])%>% # create a df of just the NLCD GAGES II predictors,
  mutate(across(where(is.numeric), ~./100))%>% # convert to percent between 0-1 to match CDL reclass
  mutate(across(where(is.numeric), round, 3))%>% #round 
  mutate(G2_Wetlands_eq = WOODYWETNLCD06+EMERGWETNLCD06,# create combinaitons for pasture and wetlands
         .keep = 'unused')%>%
  left_join(., df.NWIS.CDL.reclass%>%select(Name, Developed, Forest, Ag, Water,Pasture, Wetlands_all, Other), by = c('STAID'='Name')) %>% # merge with CDL reclass to compare (and rearrange the columns of the right joined df):
  mutate(Dev_diff = Developed - DEVNLCD06,
         Forest_diff = Forest - FORESTNLCD06,
         Ag_diff = Ag - PLANTNLCD06,
         Water_diff =  Water - WATERNLCD06,
         Pasture_diff = Pasture - PASTURENLCD06,
         Wetlands_diff = Wetlands_all - G2_Wetlands_eq,
         .keep = 'unused')

temp%>%
  pivot_longer(cols = -STAID)%>%
  ggplot(., aes(x = STAID, y = name, fill = value)) +
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "yellow")

# Even though I tried to tweak the reclassification matrix I am going to do NLCD too:






















####~~~~ Recreating the GAGES II predictor set ~~~~####





#### Land Use: NLCD ####

# the CDL and NLCD differ in that the CDL combined Pasture and grassland while the nLCD does not
# thus, even though the CDL can be aggregated to et close to the NLCD, I want to also have the nLCD
# to the difference.

# download: to do this:

# NLCD is not working when trying to download based on the entire polygon df, so going to use lapply to download individually:
# 
# l.rast.NWIS.NLCD.2006 <- lapply(seq_along(df.sf.NWIS$Name), \(i) get_nlcd(template = st_cast(df.sf.NWIS, "MULTIPOLYGON")[i,], label = as.character(i), year = 2006))
# l.rast.NWIS.NLCD.2016 <- lapply(seq_along(df.sf.NWIS$Name), \(i) get_nlcd(template = st_cast(df.sf.NWIS, "MULTIPOLYGON")[i,], label = as.character(i), year = 2016))

# convert to SpatRasters:

# l.rast.NWIS.NLCD.2006<-lapply(l.rast.NWIS.NLCD.2006, rast)
# l.rast.NWIS.NLCD.2016<-lapply(l.rast.NWIS.NLCD.2016, rast)

# save SpatRasters as tif files:

# lapply(seq_along(l.rast.NWIS.NLCD.2006), \(i) writeRaster(l.rast.NWIS.NLCD.2006[[i]], filename = paste0('Downloaded_Data/NLCD_rasters/NLCD_2006_SpatRaster_for', df.sf.NWIS$Name[i], '.tif'), overwrite=TRUE))
# lapply(seq_along(l.rast.NWIS.NLCD.2016), \(i) writeRaster(l.rast.NWIS.NLCD.2016[[i]], filename = paste0('Downloaded_Data/NLCD_rasters/NLCD_2016_SpatRaster_for', df.sf.NWIS$Name[i], '.tif'), overwrite=TRUE))

# create vectors of raster names 

names.2006<-paste0('Downloaded_Data/NLCD_rasters/NLCD_2006_SpatRaster_for', df.sf.NWIS$Name, '.tif')
names.2016<-paste0('Downloaded_Data/NLCD_rasters/NLCD_2016_SpatRaster_for', df.sf.NWIS$Name, '.tif')

# load SpatRasters back in:

# l.rast.NWIS.NLCD.2006<-lapply(names.2006, rast)
# l.rast.NWIS.NLCD.2016<-lapply(names.2016, rast)

# see if 2006 and 2016 crs are the same:

# crs(l.rast.NWIS.NLCD.2006[[2]])
# crs(l.rast.NWIS.NLCD.2016[[2]])

# yes

# reproject to sample watershed vector data to match raster data:

# vect.NWIS<-vect(df.sf.NWIS) # need to first ininalize vect since sometimes reading in rdata file (terra issue)
# 
# vect.NWIS.proj<-terra::project(vect.NWIS, crs(l.rast.NWIS.NLCD.2006[[1]]))

# extract frequency tables for each sample watershed

# system.time({l.NWIS.NLCD.2006 <- lapply(seq_along(l.rast.NWIS.NLCD.2006), \(i) terra::extract(l.rast.NWIS.NLCD.2006[[i]], vect.NWIS.proj[i], ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2)))})
# system.time({l.NWIS.NLCD.2016 <- lapply(seq_along(l.rast.NWIS.NLCD.2016), \(i) terra::extract(l.rast.NWIS.NLCD.2016[[i]], vect.NWIS.proj[i], ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2)))})

# save(l.NWIS.NLCD.2006, file = 'Processed_Data/l.NWIS.NLCD.2006.Rdata')
# save(l.NWIS.NLCD.2016, file = 'Processed_Data/l.NWIS.NLCD.2016.Rdata')

load("Processed_Data/l.NWIS.NLCD.2006.Rdata")
load("Processed_Data/l.NWIS.NLCD.2016.Rdata")

# reclassify: the GAGES II predictors for NLCD land use are the sum of a few NLCD classes (see the kable table of the variable descriptions made above). To do this:

# create and ordered vector based on legend (legend comes from Ryan_funcitons.R)

Class3.for.G2<-c("WATERNLCD06", "SNOWICENLCD06", rep("DEVNLCD06", 4), "BARRENNLCD06", "DECIDNLCD06", "EVERGRNLCD06", "MIXEDFORNLCD06", NA, "SHRUBNLCD06", "GRASSNLCD06", NA, NA, NA, "PASTURENLCD06", "CROPSNLCD06", "WOODYWETNLCD06", "EMERGWETNLCD06")

# add an identifier to this vector so the column names are slightly different than the GAGES II predictors:

Class3.for.G2<-paste0('R_', Class3.for.G2)

# create a df from this vector (for latter use):

df.Class3<-data.frame(Class = unique(Class3.for.G2)[complete.cases(unique(Class3.for.G2))])
  
# create new df from the legend dataframe:

legend.for.G2<-legend%>%mutate(Class3 = Class3.for.G2)

# reclassify the NLCD using this new legend and clean up the dataframe from the next step

l.NWIS.NLCD.2006<-lapply(l.NWIS.NLCD.2006, \(i) left_join(as.data.frame(i), legend.for.G2%>%select(Class, Class3), by = 'Class')%>%mutate(Class = Class3)%>%select(-Class3)%>%group_by(Class)%>%summarise(Freq = sum(Freq)))
l.NWIS.NLCD.2016<-lapply(l.NWIS.NLCD.2016, \(i) left_join(as.data.frame(i), legend.for.G2%>%select(Class, Class3), by = 'Class')%>%mutate(Class = Class3)%>%select(-Class3)%>%group_by(Class)%>%summarise(Freq = sum(Freq)))

# add back missing varaibles using df.Class3:

l.NWIS.NLCD.2006<-lapply(l.NWIS.NLCD.2006, \(i) left_join(df.Class3, i, by = 'Class')%>%replace(is.na(.), 0))
l.NWIS.NLCD.2016<-lapply(l.NWIS.NLCD.2016, \(i) left_join(df.Class3, i, by = 'Class')%>%replace(is.na(.), 0))

# pivot_wider the df in the lists,
# add a Name column for the site,
# and add new columns for FOREST, and PLANT:

l.NWIS.NLCD.2006<-lapply(seq_along(l.NWIS.NLCD.2006), \(i) l.NWIS.NLCD.2006[[i]]%>%
                           group_by(Class)%>%
                           summarise(Freq = sum(Freq))%>%
                           pivot_wider(names_from = Class, values_from = Freq)%>%
                           mutate(Name = df.sf.NWIS$Name[i], .before = 1)%>%
                           mutate(R_FORESTNLCD06 = R_DECIDNLCD06+R_EVERGRNLCD06+R_MIXEDFORNLCD06,
                                  R_PLANTNLCD06 = R_PASTURENLCD06+R_CROPSNLCD06)%>%
                           as.data.frame(.))
l.NWIS.NLCD.2016<-lapply(seq_along(l.NWIS.NLCD.2016), \(i) l.NWIS.NLCD.2016[[i]]%>%
                           group_by(Class)%>%
                           summarise(Freq = sum(Freq))%>%
                           pivot_wider(names_from = Class, values_from = Freq)%>%
                           mutate(Name = df.sf.NWIS$Name[i], .before = 1)%>%
                           mutate(R_FORESTNLCD06 = R_DECIDNLCD06+R_EVERGRNLCD06+R_MIXEDFORNLCD06,
                                  R_PLANTNLCD06 = R_PASTURENLCD06+R_CROPSNLCD06)%>%
                           as.data.frame(.))

# bind the lists into a single dataframe
# note some of the sites have different length dtaframes because they didnt have all the same number of NLCD classes. When binding rows this will give a dataframe of the maximum length and put NAs for sites where there wasn't a column: 

df.NWIS.NLCD.2006<-bind_rows(l.NWIS.NLCD.2006)
df.NWIS.NLCD.2016<-bind_rows(l.NWIS.NLCD.2016)

# now compare to GAGES II:

df.compare.G2to2006<-df.G2.reduced%>%select(STAID, names(df.G2.reduced)[c(8:11,18,19,21,22)])%>% # create a df of just the NLCD GAGES II predictors,
  mutate(across(where(is.numeric), ~./100))%>%
  mutate(across(where(is.numeric), round, 3))%>%
  # mutate_all(~ ifelse(is.numeric(.), round(.*0.01, 2), .))%>%
  left_join(., df.NWIS.NLCD.2006, by = c('STAID'='Name'))%>%
  pivot_longer(cols = -STAID)%>%
  mutate(name = sub('*R_', '', name))%>%
  group_by(STAID, name) %>%
  summarise(diff = diff(value)) %>%
  pivot_wider(names_from = name, values_from = diff) %>%
  rename_at(-1, ~paste0(., "_diff"))%>%
  mutate(across(where(is.numeric), round, 3))

df.compare.G2to2016<-df.G2.reduced%>%select(STAID, names(df.G2.reduced)[c(8:11,18,19,21,22)])%>% # create a df of just the NLCD GAGES II predictors,
  mutate(across(where(is.numeric), ~./100))%>%
  mutate(across(where(is.numeric), round, 3))%>%
  # mutate_all(~ ifelse(is.numeric(.), round(.*0.01, 2), .))%>%
  left_join(., df.NWIS.NLCD.2016, by = c('STAID'='Name'))%>%
  pivot_longer(cols = -STAID)%>%
  mutate(name = sub('*R_', '', name))%>%
  group_by(STAID, name) %>%
  summarise(diff = diff(value)) %>%
  pivot_wider(names_from = name, values_from = diff) %>%
  rename_at(-1, ~paste0(., "_diff"))%>%
  mutate(across(where(is.numeric), round, 3))

# plot as heatmap:

df.compare.G2to2006%>%
  pivot_longer(cols = -STAID)%>%
  ggplot(., aes(x = STAID, y = name, fill = value)) +
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "yellow")

df.compare.G2to2016%>%
  pivot_longer(cols = -STAID)%>%
  ggplot(., aes(x = STAID, y = name, fill = value)) +
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "yellow")

#### Land Use: CDL ####

# the CDL ones in GAGES II do not match the CDL names :/
# to match them up:

# get the names of the crop from the GAGESII CDL ones:

v<-names(df.G2.reduced)
cdl_names <- grep("^CDL", v)
v <- v[cdl_names]
v_suffix<-sub("^CDL_", "", v)

# now look at the two and compare:

v_suffix # GAGES II suffix

names(df.NWIS.CDL.2008) # CDL datalayers

# Soybeans, Corn, Alfalfa, Spring_Wheat, Other_Hay/Non_Alfalfa, Grassland/Pasture, Dbl_Crop_WinWht/Soybeans, Winter_Wheat
# *note*: ALL_OTHER_LAND is non-crop land, so excluding

keep<-c('Name', 'Corn', 'Soybeans', 'Spring_Wheat', 'Winter_Wheat', 'Dbl_Crop_WinWht/Soybeans', 'Alfalfa', 'Other_Hay/Non_Alfalfa', 'Grassland/Pasture')

# now extract just these from the CDL datalayers df:

temp<-df.NWIS.CDL.2008%>%select(keep)

# rename them to the GAGES II CDL names with a unique ID for it being the datalayers one ([-9] to remove ALL_OTHER_LAND):

names(temp)<-paste0('R_', c('Name', v[-9]))

# now add these columns to the GAGES II df to compare:
# also set NA to zero and round to 2 decimial places:

df.compare.G2toCDL<-left_join(df.G2.reduced%>%select(STAID, v[-9])%>%mutate(across(where(is.numeric), ~./100)), temp, by = c('STAID'='R_Name'))%>%
  replace(is.na(.), 0)%>%
  mutate(across(where(is.numeric), round, 2))

# now compare GAGES II to CDL:

df.compare.G2toCDL<-df.compare.G2toCDL%>%
  pivot_longer(cols = -STAID)%>%
  mutate(name = sub('*R_', '', name))%>%
  group_by(STAID, name) %>%
  summarise(diff = diff(value)) %>%
  pivot_wider(names_from = name, values_from = diff) %>%
  rename_at(-1, ~paste0(., "_diff"))

# and heat map:

df.compare.G2toCDL%>%
  pivot_longer(cols = -STAID)%>%
  ggplot(., aes(x = STAID, y = name, fill = value)) +
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "yellow")

#




















#### Elevation: NED ####

# download:

# DEM.NWIS<-get_ned(df.sf.NWIS, label = '2') # already SpatRaster!

# writeRaster(DEM.NWIS, file='Downloaded_Data/DEM.NWIS.tif', overwrite=TRUE)

# load in rast:

DEM.NWIS<-rast('Downloaded_Data/DEM.NWIS.tif')

# plot(DEM.NWIS)

# extract elevation metrics over each sample watershed: to do this:
# build a function with multiple functions:

# f <- function(x, na.rm = T) {
#   c(mean=mean(x, na.rm = na.rm),
#     median=median(x, na.rm = na.rm),
#     range=max(x, na.rm = na.rm)-min(x, na.rm = na.rm),
#     sd=sd(x, na.rm = na.rm)
#   )
# }

# reproject NWIS basins to DEM crs:

# vect.NWIS<-vect(df.sf.NWIS)

# vect.NWIS.proj<-terra::project(vect.NWIS, crs(DEM.NWIS))

# extract the metrics over each watershed using the function above:

# df.NWIS.DEM <- as.data.frame(terra::extract(DEM.NWIS, vect.NWIS.proj, f))

# set the names of the df:

# names(df.NWIS.DEM)<-c('Name', 'Elev_Avg','Elev_Median', 'Elev_Range', 'Elev_SD')

# set the names of the sites:

# df.NWIS.DEM$Name<-df.sf.NWIS$Name

# finally save the df:

# save(df.NWIS.DEM, file= 'Processed_Data/df.NWIS.DEM.Rdata')

load('Processed_Data/df.NWIS.DEM.Rdata')

# compare to GAGES II:

# Only the median and STD are in GAGES II:

df.compare.G2toNED<-left_join(df.NWIS.DEM%>%select(c(1,3,5))%>%rename(R_ELEV_MEDIAN_M_BASIN = 2, R_ELEV_STD_M_BASIN = 3), df.G2.reduced%>%select(STAID, ELEV_MEDIAN_M_BASIN,ELEV_STD_M_BASIN), by = c('Name'='STAID'))%>%
  mutate(across(where(is.numeric), round, 0))%>%
  pivot_longer(cols = -Name)%>%
  mutate(name = sub('*R_', '', name))%>%
  group_by(Name, name) %>%
  summarise(diff = diff(value)) %>%
  pivot_wider(names_from = name, values_from = diff) %>%
  rename_at(-1, ~paste0(., "_diff"))
  
# make heat map:

df.compare.G2toNED%>%
  pivot_longer(cols = -Name)%>%
  ggplot(., aes(x = Name, y = name, fill = value)) +
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "yellow")

# done with DEM





























#### Soils: SURRGO ####

# example workflow:

# created function to download and get watershed percent of HSG:
# the funciton takes a single site from a sf.df, so use it in lapply for all sites:
# want to test two random sites to GAGES 2 and streamstats HSG values:

# get two random sites from df.sf.NWIS:

keep<-sample(df.sf.NWIS$Name, 2)

# check map (to make sure the sites arent too big since I need to download soils data for them):

x<-df.sf.NWIS%>%filter(Name %in% keep)
mapview(x, zcol = 'Name')

# looks like two good random sites

# get SURRGO for these sites:

l.test.soils<-lapply(keep, \(i) fun.SURRGO_HSG(i, df.sf.NWIS))

names(l.test.soils)<-keep

# look at how much of the watershed is represented:

lapply(l.test.soils, \(i) sum(i$Watershed_Percent))

df.test.soils<-bind_rows(l.test.soils, .id = 'STAID')%>%filter(HSG %in% c('A', 'B', 'C', 'D'))%>%rename(Datalayers = 3)

# QA with GAGES II:

x<-df.G2.reduced%>%
  select(STAID, c('HGA', 'HGB', 'HGC', 'HGD'))%>%
  mutate(across(where(is.numeric), ~./100))%>%
  rename('A' = 2, 'B'=3, 'C'=4, 'D'=5)%>%
  filter(STAID %in% keep)%>%
  pivot_longer(cols=-STAID, names_to = 'HSG', values_to = 'G2')%>%
  left_join(., df.test.soils, by = c('STAID', 'HSG'))%>%
  mutate(Diff = Datalayers-G2)

# QA that dfwith streamstats (in df.sf.NWIS):

y<-df.sf.NWIS%>%filter(Name %in% keep)%>%select(Name, starts_with('SSURGO'))

# discrepancy!

# full workflow:

# get SURRGO for these sites:

l.soils<-lapply(df.sf.NWIS$Name, \(i) fun.SURRGO_HSG(i, df.sf.NWIS))

names(l.test.soils)<-keep

# look at how much of the watershed is represented:

lapply(l.test.soils, \(i) sum(i$Watershed_Percent))

df.test.soils<-bind_rows(l.test.soils, .id = 'STAID')%>%filter(HSG %in% c('A', 'B', 'C', 'D'))%>%rename(Datalayers = 3)

























#### Roads ####

# install.packages("tigris")

# library(tigris)

# download roads:

# NY.roads<-roads(state = 36, county = 'Onondaga')

# mapview(NY.roads)

# I dont think it is worth processing the roads right now...



























#### NHD ####

# install.packages("nhdplusTools")

library(nhdplusTools)
library(sf)

# read in NWIS site outlet corrdinate:

# df.NWIS.TP_site_metadata<-read.csv("Raw_Data/df.NWIS.TP_site_metadata.csv", colClasses = c(site_no = "character"))%>%
#   filter(site_no %in% df.sf.NWIS$Name)

# order df to match that of names.2006 (from df.sf.NWIS) for use in function below:

# df.site_metadata<-df.NWIS.TP_site_metadata[match(df.sf.NWIS$Name, df.NWIS.TP_site_metadata$site_no),]

# use function to get 100 and 800 meter rip buffer percent dev, plant, and forest:

# l.RIP<-lapply(seq_along(df.sf.NWIS$Name), \(i) fun.NHD.RIP.buffers(i, df.site_metadata))

# combine list into single df:

# df.RIP<-bind_rows(l.RIP)%>%mutate(Name = df.sf.NWIS$Name)

# save df:

# save(df.RIP, file='Processed_Data/df.NHD.RIP.buffer.100.800.pland.Rdata')

load('Processed_Data/df.NHD.RIP.buffer.100.800.pland.Rdata')

# QA with gauges 2: these should line up with the RIP variables:
  
df.compare.G2toNHD<-df.G2.reduced%>%
  select(STAID, contains('RIP'))%>%
  mutate(across(where(is.numeric), ~./100))%>%
  mutate(across(where(is.numeric), round, 3))%>%
  left_join(., df.RIP, by = c('STAID'='Name'))%>%
  pivot_longer(cols = -STAID)%>%
  mutate(name = sub('*R_', '', name))%>%
  group_by(STAID, name) %>%
  summarise(diff = diff(value)) %>%
  pivot_wider(names_from = name, values_from = diff) %>%
  rename_at(-1, ~paste0(., "_diff"))

df.compare.G2toNHD%>%
  pivot_longer(cols = -STAID)%>%
  ggplot(., aes(x = STAID, y = name, fill = value)) +
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "yellow")
  
# wow, they check out pretty good!

#

####~~~~Finalize DataLayers ~~~~####


df.NWIS.NLCD.2016

df.NWIS.CDL.2008

df.NWIS.DEM





# map:

library(maptiles)
bg <- get_tiles(ext(vect.flowline.buffer.100))

plotRGB(bg)
lines(vect.flowline.buffer.100, col="blue", lwd=3)





























