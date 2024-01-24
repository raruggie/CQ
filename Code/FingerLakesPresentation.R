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

# 1)
# code of figures and any other anaylsis for fingers lakes presentation

####################### Workflow #######################

#### Map of sites ####

# read in df.sf.NWIS:

load('Processed_Data/NWIS_Watershed_Shapefiles.Rdata')

# arrange df by drainage area sizeto help with plotting order in mapview:

df.sf.NWIS<-arrange(df.sf.NWIS, desc(drain_area_va))

# download metadata for these sites:

df.points<-readNWISsite(df.sf.NWIS$Name)%>%
  st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4326)%>%
  arrange(desc(drain_area_va))

# make map:

mapview(df.sf.NWIS, zcol = 'drain_area_va', layer.name = 'Drainage Area')+mapview(df.points, zcol = 'drain_area_va', legend = FALSE)

#

#### CQ curves ####

# load in plotting df:

load('Processed_Data/TP.df_Seg.2.Rdata')

# make map of sites colored by land use: to do this first need to add the land use class to df.sf.NWIS:

df.sf.NWIS<-df.sf.NWIS%>%left_join(., df_Seg.2%>%distinct(site, .keep_all = T)%>%select(site, USGS.LU.Adjusted, Type), by = c('Name'='site'))
df.points<-df.points%>%left_join(., df_Seg.2%>%distinct(site, .keep_all = T)%>%select(site, USGS.LU.Adjusted, Type), by = c('site_no'='site'))

mapview(df.sf.NWIS, zcol = 'USGS.LU.Adjusted', layer.name = 'USGS Landuse Classification')+mapview(df.points, zcol = 'drain_area_va', legend = FALSE)

# create 3 pairs of sites:

keep<-c("04249000", "01357500", "01362497", "04231000", "04269000", "04232050")
keep_names<-c("Oswego", 'Mohawk', 'LBK', 'Black', 'SRR', 'Allen')
df_keep<-data.frame(keep, keep_names)

# filter df.Seg_2 to 3 pairs of sites:

df_Seg.2<-filter(df_Seg.2, site %in% keep)

# and add site name columns:
# first change the string for little beaver kill:

df.points$station_nm[31]<-"LITTLE BEAVER KILL AT\nBEECHFORD NEAR MT TREMPER NY"

df_Seg.2<-left_join(df_Seg.2, df.points%>%select(site_no,station_nm, drain_area_va), by = c('site'='site_no'))%>%mutate(site = paste(site, station_nm, '\nDA =', drain_area_va, 'sqmi'))

# set site column as ordered factor based on DA size for facet plotting order:

df_Seg.2<-df_Seg.2%>%mutate(site = factor(site, levels = unique(site[order(drain_area_va)])))

# try plotting all together:

p<-ggplot(df_Seg.2, aes(x = log(Q_real), y = log(C)))+
  geom_point(aes(color = Type))+
  scale_color_manual(name = "CQ Type", values = c("red", "blue", "green"))+
  ylab('log(TP - mg/L)')+
  xlab('log(Discharge - cfs)')+
  geom_smooth(method = 'lm')+
  new_scale_color() +
  # geom_line(aes(x = Q, y = Seg_C), size = 2.5, color = 'black')+
  # geom_line(aes(x = Q, y = Seg_C, color = slope_angle), size = 2)+
  # scale_color_manual(name = "Slope Angle", values = hc)+
  facet_wrap('site', scales = 'fixed', ncol = 2, dir="v")+
  theme(
    # strip.background = element_blank(),
    # strip.text.x = element_blank(),
    legend.position="bottom"
  )+
  geom_rect(data = df_Seg.2%>%distinct(df_Seg.2$site, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = USGS.LU.Adjusted), alpha = .35)+
  scale_fill_manual(name = "USGS Landuse\n(Adjusted)", values = c("red", "blue","yellow", "green"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

p

#### Traingle plots ####

# read in df.tri for TP, TN, TDP, and SRP and rename:

load('Processed_Data/df.tri.Rdata')
df.tri.TP<-df.tri
load('Processed_Data/df.tri.SRP.Rdata')
df.tri.SRP<-df.tri
load('Processed_Data/df.tri.TN.Rdata')
df.tri.TN<-df.tri
load('Processed_Data/df.tri.TDP.Rdata')
df.tri.TDP<-df.tri

# merge site name:

df.tri.TP<-left_join(df.tri.TP, df.points%>%select(site_no,station_nm, drain_area_va), by = 'site_no')%>%mutate(site = paste(site_no, station_nm, '\nDA =', drain_area_va, 'sqmi'))
df.tri.SRP<-left_join(df.tri.SRP, df.points%>%select(site_no,station_nm, drain_area_va), by = 'site_no')%>%mutate(site = paste(site_no, station_nm, '\nDA =', drain_area_va, 'sqmi'))
df.tri.TDP<-left_join(df.tri.TDP, df.points%>%select(site_no,station_nm, drain_area_va), by = 'site_no')%>%mutate(site = paste(site_no, station_nm, '\nDA =', drain_area_va, 'sqmi'))
df.tri.TN<-left_join(df.tri.TN, df.points%>%select(site_no,station_nm, drain_area_va), by = 'site_no')%>%mutate(site = paste(site_no, station_nm, '\nDA =', drain_area_va, 'sqmi'))

# create facet plot for all 4: to do this:

# create a list of the dfs:

l.tri<-list(df.tri.TP, df.tri.TN, df.tri.TDP, df.tri.SRP)%>%purrr::set_names(c('TP', 'TN', 'TDP', 'SRP'))

# scale each slope magintude for each consituent: to do this:
# first take abs, then normalize using function in Code/Ryan_functions.R:

l.tri<-lapply(l.tri, \(i) i%>%mutate(Slope = normalized(abs(Slope))))

# bind into single df and add column with consituent id (name of list element works here):
# also highlight a few sites:

df.tri.4<-bind_rows(l.tri, .id = 'Consit')%>%
  mutate(group = ifelse(site_no %in% c("04231600","04260500","04249000","01357500"), site_no, NA))%>%
  left_join(., df_keep, by = c('site_no'='keep'))
# ready to create facet plot:

p.4<-ggplot(df.tri.4, aes(x=DEVNLCD06, y=PLANTNLCD06*10
                            , label = keep_names
                            )) +
  # geom_point(data=df.tri.4[df.tri.4$group == "important",],color="red",size=5)+
  geom_point(aes(shape=as.factor(Type), fill=abs(Slope)), size = 4) +
  scale_shape_manual(values = c("Mobilization" = 24, "Dilution" = 25, "Stationary" = 22),
                     guide = guide_legend(override.aes = list(fill = "pink")))+
  scale_fill_gradient(low = "yellow", high = "red")+
  facet_wrap('Consit', scales = 'fixed')+
  theme(
    # strip.background = element_blank(),
    # strip.text.x = element_blank(),
    legend.position="bottom"
  )+
  ggrepel::geom_text_repel(position = position_nudge(x = .7, y = .7), angle = 45, hjust = 0, vjust = 0)+
  xlab('Percent Developed')+
  ylab('Percent Agriculture')+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_y_continuous(expand = expansion(mult = .2)) # +
  # ggtitle(paste('TP', 'CQ type and slope magnitude as a function of percent Ag and Developed Land'))

p.4$labels$fill <- "Normalized (0-1) \nSlope Magnitude"
p.4$labels$shape <- "CQ Type"

p.4

# look at sites that are dilutionary:

df.tri.TP$site_no[df.tri.TP$Type=='Dilution'] # there is only 1 TP site

df.tri.SRP$site_no[df.tri.SRP$Type=='Dilution'] # there are 7 SRP sites

# make a map of the dilutionary sites:

x<-df.sf.NWIS%>%filter(Name %in% df.tri.SRP$site_no[df.tri.SRP$Type=='Dilution'])
y<-df.points%>%filter(site_no %in% df.tri.SRP$site_no[df.tri.SRP$Type=='Dilution'])

mapview(x, zcol = 'USGS.LU.Adjusted', layer.name = 'Drainage Area') #+mapview(df.points, zcol = 'drain_area_va', legend = FALSE)

# make a map of all sites with CQ trend for SRP: to do this: need to create new dfs for df.sf.NWIS and df.points with Type column for SRP data (at this point in the code thesedataframes are using TP Type):

x<-df.sf.NWIS%>%select(-Type)%>%left_join(., df.tri.SRP%>%select(site_no, Slope, Type), by = c('Name'='site_no'))%>%drop_na(Type)
y<-df.points%>%select(-Type)%>%left_join(., df.tri.SRP%>%select(site_no, Slope, Type), by = 'site_no')%>%drop_na(Type)

mapview(x, zcol = 'Type', layer.name = 'SRP Export Regime')+mapview(y, zcol = 'Type', legend = FALSE)

# what are the dilutionary SRP sites doing for TP?:

x<-df.tri.TP%>%filter(site_no %in% df.tri.SRP$site_no[df.tri.SRP$Type=='Dilution'])

#### Combining MLR ####

# the tables comparing the mLR models are for thedifferent CQ parameters for a single consitient
# for the poster I want a table of the same parameter for different consiuents:

# read in the m.lists from NWIS_X.R files and rename:

load('Processed_Data/m.list.TP.Rdata')
m.list.TP<-m.list
load('Processed_Data/m.list.TN.Rdata')
m.list.TN<-m.list
load('Processed_Data/m.list.TDP.Rdata')
m.list.TDP<-m.list
load('Processed_Data/m.list.SRP.Rdata')
m.list.SRP<-m.list

# make a list of these lists:

l.m.list<-list(m.list.TP, m.list.TN, m.list.TDP, m.list.SRP)%>%purrr::set_names(c('TP', 'TN', 'TDP', 'SRP'))

# extract the first and last element of the list to get a list of the OLS intercept and yield models for each constituent respectively:

l.m.list.intercept<-lapply(l.m.list, \(i) i[[1]])
l.m.list.yield<-lapply(l.m.list, \(i) i[[5]])

# now use tab_model to compare and make nice tables!:

tab_model(l.m.list.intercept, dv.labels = names(l.m.list.intercept), title = paste('Comparison of MLR models for CQ Intercept'), file="temp.html")
tab_model(l.m.list.yield, dv.labels = names(l.m.list.yield), title = paste('Comparison of MLR models for AAY'), file="temp.html")


#### Kable table of site result list ####

# read in dataframe of site result list:

x<-read.csv("Processed_Data/NWIS_query_result_table_for_poster.csv", check.names=FALSE)

x %>%
  kbl(align = "c",escape = F, caption = "HIIIII") %>%
  kable_classic(html_font = 'Times', font_size = 14, full_width = F) %>%
  add_header_above(c(" ", "Number of Sites"=4)) 


%>%
  add_header_above(c(" " = 2, "2019" = 4, " " = 1, "2020" = 4, " " = 1, "2021" = 4, " " = 1)) %>%
  # add_header_above(c(" " = 1, "TP Loading (g/ha) by Water Year (Oct-Sep)\nG = Growing Season (May-Sep), NG = Non-growing" = 15))%>%
  collapse_rows(., columns = 1, valign = 'middle')%>%
  row_spec(c(3,6), extra_css = "border-bottom: 1px solid")%>%
  kableExtra::footnote(general = "Sampling errors resulted in the following number of missing events and load estimates:", 
                       alphabet = c(
                         "Growing: 1 event, 171 g/ha, Non-Growing: 4 events, 261 g/ha",
                         "Growing: 3 events, 64 g/ha, Non-Growing: 7 events, 630 g/ha",
                         "Non-Growing: 14 events, 2770 g/ha",
                         "Growing: 1 event, 2 g/ha",
                         "Non-Growing: 14 events, 555 g/ha",
                         "Growing: 1 event, 3 g/ha"
                       )
  )

# 














