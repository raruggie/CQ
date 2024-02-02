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

# load in consitieunt CQ dfs

load("Processed_Data/TP.Rdata")
load("Processed_Data/TN.Rdata")
load("Processed_Data/NO3.Rdata")
load("Processed_Data/TDP.Rdata")
load("Processed_Data/SRP.Rdata")

# create list of these dfs:

l<-grep("df.",names(.GlobalEnv),value=TRUE)
l<-do.call("list",mget(l))
l<-Filter(is.list, l)

# extract just the unique site numbers from each df:

l<-lapply(l, \(i) unique(i$site_no))

# perform intersection across list:

Reduce(intersect, l) # zero!!

# # TP, TN, SRP, TDP (i.e. removing nitrate)

l['df.NO3_CQ']<-NULL

Reduce(intersect, l) # 4

# overlap between some combinations:

# TP and TN:

Reduce(intersect, l[c("df.TP_CQ", "df.TN_CQ")]) # 17

# TP and SRP:

Reduce(intersect, l[c("df.TP_CQ", "df.SRP_CQ")]) # 40

# TP and TDP:

Reduce(intersect, l[c("df.TP_CQ", "df.TDP_CQ")]) # 18

# TP, TN, and SRP:

Reduce(intersect, l[c("df.TP_CQ", "df.TN_CQ", "df.SRP_CQ")]) # also 17

# Looking at these 17 sites:
# from each NWIS .R code, I need to export the df_Seg.2 dataframe and read in here:

load('Processed_Data/TP.df_Seg.2.Rdata')
TP.df_Seg.2<-df_Seg.2
load('Processed_Data/TN.df_Seg.2.Rdata')
TN.df_Seg.2<-df_Seg.2
load('Processed_Data/SRP.df_Seg.2.Rdata')
SRP.df_Seg.2<-df_Seg.2

# make list of these dataframes:

l.df_Seg.2<-list(TP.df_Seg.2, TN.df_Seg.2, SRP.df_Seg.2)%>%purrr::set_names(c('TP', 'TN', 'SRP'))

# filter each dataframe to just the 17 sites:

l.df_Seg.2<-lapply(l.df_Seg.2, \(i) i%>%filter(site %in% Reduce(intersect, l[c("df.TP_CQ", "df.TN_CQ", "df.SRP_CQ")])))

# the n_sample_rank column is used in the facet plot function for ordering the sites
# but between consitents they might have different n_sample_ranks
# thus I am going to set the site column as an ordered factor for the TP dataframe,
# then use that to set the other two consitents site orders so they align

# Convert the character column to an ordered factor based on numeric_column

l.df_Seg.2[['TP']]$site <- factor(l.df_Seg.2[['TP']]$site, levels = unique(l.df_Seg.2[['TP']]$site)[order(l.df_Seg.2[['TP']]$n_sample_rank)], ordered = TRUE)

# check to see if this worked:

y<-distinct(l.df_Seg.2[['TP']], site, .keep_all = T) # viewing this dataframe, thesitecolumn is in the order of the n_sample_rank column

y$site==levels(l.df_Seg.2[['TP']]$site) # looks good

# Use the levels from TP df to set the levels of TN and SRP dfs:

l.df_Seg.2[['TN']]$site <- factor(l.df_Seg.2[['TN']]$site, levels = levels(l.df_Seg.2[['TP']]$site), ordered = TRUE)

l.df_Seg.2[['SRP']]$site <- factor(l.df_Seg.2[['SRP']]$site, levels = levels(l.df_Seg.2[['TP']]$site), ordered = TRUE)

# the Type column is not represented fully across all constituents
# It is only fully represented in the SRP df.
# so doing a similar thing with the type column

l.df_Seg.2[['SRP']]$Type <- factor(l.df_Seg.2[['SRP']]$Type, levels = unique(l.df_Seg.2[['SRP']]$Type), ordered = TRUE)
l.df_Seg.2[['TP']]$Type <- factor(l.df_Seg.2[['TP']]$Type, levels = levels(l.df_Seg.2[['SRP']]$Type), ordered = TRUE)
l.df_Seg.2[['TN']]$Type <- factor(l.df_Seg.2[['TN']]$Type, levels = levels(l.df_Seg.2[['SRP']]$Type), ordered = TRUE)

# now make plots for each:

l.plots<-list()

for (i in names(l.df_Seg.2)){
  
  hc<-heat.colors(length(unique(l.df_Seg.2[[i]]$slope_angle)), rev = T)
  
  p<-ggplot(l.df_Seg.2[[i]], aes(x = log(Q_real), y = log(C)))+
    geom_point(aes(color = Type))+
    scale_color_manual(name = "CQ Type", values = c("blue","green", "red"), drop = FALSE)+
    geom_smooth(method = 'lm')+
    ylab(paste0('log(', i, ')'))+
    new_scale_color() +
    geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
    geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
    scale_color_manual(name = "Slope Angle", values = hc)+
    facet_wrap(dplyr::vars(site), scales = 'free')+
    # theme(
    #   strip.background = element_blank(),
    #   strip.text.x = element_blank()
    # )+
    geom_rect(data = l.df_Seg.2[[i]]%>%distinct(l.df_Seg.2[[i]]$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .35)+
    scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","yellow", "green"))
  
  l.plots[[i]]<-p
  
}

l.plots

# the third plot for SRP looks like the TP blob prior to
# removing the pre-2001 samples, so I want to take a closer look:

x<-SRP.df_Seg.2%>%filter(site == '0422026250') # the dates look fine

# I am also thinking that this approach does not account for samples taken on the same day  
# across the sites. 

# I will use the sample dates of the TP df to filter the TN and SRP samples. TO do this:

# create a temporaty dataframe to filter with:

x<-l.df_Seg.2[['TP']]%>%mutate(x=paste(site, Date))

# filter the TN and SRP dfs but fby first creating a column to match the one in x (above):

y<-l.df_Seg.2[['TN']]%>%mutate(y=paste(site, Date))%>%filter(y %in% x$x)%>%select(-y)
z<-l.df_Seg.2[['SRP']]%>%mutate(y=paste(site, Date))%>%filter(y %in% x$x)%>%select(-y)

# compare the number of samples to orginal:

dim(y)[1]; dim(l.df_Seg.2[['TN']])[1] # three less

dim(z)[1]; dim(l.df_Seg.2[['SRP']])[1] # 124 less




















