library(sf)
library(tidyverse)
library(mapview)
library(terra)

x<-st_read('Raw_Data/Gauges_2_shapefiles/bas_nonref_NorthEast.shp')

load('Processed_Data/TP.Rdata')

x<-x%>%filter(GAGE_ID %in% df.TP_CQ$site_no)

v<-vect(x)

class(v)

plot(v)

mapview(v)
