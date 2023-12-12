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
l[[4]]<-NULL

# extract just the unique site numbers from each df:

l<-sapply(l, \(i) unique(i$site_no))

# perform intersection across list:

Reduce(intersect, l) # only get 3

# try removing nitrate:

l[[4]]<-NULL

Reduce(intersect, l) # 14

# overlap between some combinations:

# TP and TN:

Reduce(intersect, l[c(2,3)]) # 24

# TP and SRP:

Reduce(intersect, l[c(3,4)]) # 40

# TP and TDP:

Reduce(intersect, l[c(1,3)]) # 40














