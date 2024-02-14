# function to turn streamstats polygons into vects
# use 100 for small_artefacts_threshold

SS_sf_to_vect<-function(sf_list, small_artefacts_threshold){
  
  # check validity
  # lapply(sf_list, st_is_valid, reason = TRUE) |> str()
  
  # make valid, check resulting classes
  l2 <- lapply(sf_list, st_make_valid)
  # lapply(l2, st_geometry) |> 
  #   lapply(class) |>
  #   str()
  
  # some shapes ended up as multipolygons, 
  # cast all to polygons and drop small artefacts (100m^2 area threshold)
  l3 <- lapply(l2, st_geometry) |>
    lapply(st_cast, "POLYGON") |> 
    lapply(\(g) g[units::drop_units(st_area(g)) > small_artefacts_threshold])
  # lapply(l3, class) |> str()
  
  # convert to SpatVector, check validity
  l_vect <- lapply(l3, vect)
  # lapply(l_vect, is.valid) |> str()
  
  # lets deal with that as well 
  l_vect <- lapply(l_vect, makeValid)
  
  # check resulting geomtypes
  # lapply(l_vect, geomtype) |> str()
  
}

# not a function but variables used for Reassigning NLCD legend to less classes;

legend<-pal_nlcd()%>%
  mutate(ID2=c(1,2,3,3,3,3,2,4,4,4,4,4,5,5,2,2,6,7,8,8))

legend_2<-data.frame(ID2 = unique(legend$ID2))%>%
  mutate(Class2 = c("Open Water", "Other", "Developed", "Forest", "Grassland", "Pasture/Hay", "Cultivated Crops", "Wetlands"),
         ID3 = c(1,2,3,4,5,6,6,7),
         Class3 = c("Open Water", "Other", "Developed", "Forest", "Grassland", "Ag", "Ag", "Wetlands"))

legend<-left_join(legend, legend_2, by = 'ID2')

# function to determine season from date:

getSeason <- function(DATES) {
  WS <- as.Date("2012-12-15", format = "%Y-%m-%d") # Winter Solstice
  SE <- as.Date("2012-3-15",  format = "%Y-%m-%d") # Spring Equinox
  SS <- as.Date("2012-6-15",  format = "%Y-%m-%d") # Summer Solstice
  FE <- as.Date("2012-9-15",  format = "%Y-%m-%d") # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Fall")))
}

# value for scaleing ggplot colors:
cc <- scales::seq_gradient_pal("purple", "yellow", "Lab")(seq(0,1,length.out=9))

# function for elapsed_months

elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

Delineate <- function(long, lat) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      
      message("Delineating Watershed:")
      
      delineateWatershed(xlocation = long, ylocation = lat, crs = 4326, includeparameters = 'true') 
      # The return value of `streamstats::delineateWatershed` is the actual value 
      # that will be returned in case there is no condition (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message("If no other red text than it worked!")
    }
  )    
  return(out)
}

fun.l.SS_WS.to.sfdf<-function(list_of_SS_WS){
  
  list_of_SS_WS<-list_of_SS_WS[sapply(list_of_SS_WS, function(x) class(x) == "watershed")]
  
  l.sp<-list_of_SS_WS
  
  for (i in seq_along(l.sp)){
    
    tryCatch({
      
      l.sp[[i]]<-toSp(watershed = l.sp[[i]], what = 'boundary')
      
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  # set the names of the list:
  
  names(l.sp)<-names(list_of_SS_WS)
  
  # need to remove the drainage areas that did not work in toSp (nothing we can do about losing these, idk why some dont come out of the function right):
  
  l.sp<-l.sp[!sapply(l.sp, function(x) class(x) == "watershed")]
  
  # convert from sp (old gis in r) to sf (new gis in r)
  
  l.sf<-lapply(l.sp, st_as_sf)
  
  # make valid:
  
  l.sf<-lapply(l.sf, st_make_valid)
  
  # need to convert the Shape_Area column to numeric for all dfsin the list or bind_rows wont work:
  
  l.sf<-lapply(l.sf, \(i) i%>%mutate(Shape_Area = as.numeric(Shape_Area)))
  
  # create a single sf df with all the sample site draiange areas:
  
  df.sf<-bind_rows(l.sf, .id = 'Name')%>%
    relocate(Name, .before= 1)
  
  return(df.sf)
}

Ryan_toSp <- function(DA) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      
      message("converting to sp:")
      
      toSp(DA) 
      # The return value of `streamstats::delineateWatershed` is the actual value 
      # that will be returned in case there is no condition (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message("If no other red text than it worked!")
    }
  )    
  return(out)
}


# Snap a point to the closest point on a line segment using sf
# https://stackoverflow.com/questions/51292952/snap-a-point-to-the-closest-point-on-a-line-segment-using-sf

st_snap_points = function(x, y, max_dist = 1000) {
  
  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)
  
  out = do.call(c,
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  return(out)
}

# a function is created to pull the data for each consituent this will be used in lapply to create a list with dataframes for  each consituent
# whatNWISdata is again used but this time the pcode for the different nutrients are inputs to the function created below:

# test variables for function building:

# parameter_code<-'00600' # pcode for TP
# 
# df.flow_query<-df.NWIS.Q_sites
# 
# state<-'NY'

fun.df.Pair_consit_flow<-function(parameter_code, df.flow_query, n_samples = 20, state = 'NY') {
  
  # download metadata for the sites in NY with data of the constituent of interest:
  
  df.consit_fun<-whatNWISdata(stateCd=state,parameterCd=parameter_code)
  
  # reduce duplicate site entries:
  
  df.consit_fun<-df.consit_fun%>% 
    group_by(site_no)%>% 
    slice(which.max(count_nu)) 
  
  # make sure the site_no columns are the same class:
  
  df.flow_query<-df.flow_query%>%mutate(site_no = as.character(site_no))
  
  # pair up this dataframe with the flow query dataframe:
  # select only asubsetof columns
  # filter to the sites returned in the flow query
  # then left join with the flow query df (could have just done a right join instead of filter+left_join)
  # filter for sites with over 99 samples
  # arrange the resulting df by the number of samples and the number of flow observations:
  
  df.consit_flow_fun<-df.consit_fun%>%
    select(c(2,3,5,6,22,23,24))%>% # "site_no"     "station_nm"  "dec_lat_va"  "dec_long_va" "begin_date"  "end_date"    "count_nu"
    filter(site_no %in% df.flow_query$site_no)%>%
    left_join(.,df.flow_query[,c(2,3,5,6,22,23,24)])%>%
    filter(count_nu > n_samples)%>%
    arrange(desc(count_nu), desc(nflowdays))
  
  # return this dataframe as the result of this function
  return(df.consit_flow_fun)
  

}

# read in allsheets of a excel workbook into a list (used it for gauges2)

library(readxl)    

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# function to take dataframe(s) of 21 day climate lags and turn into a single 1 row df of the processed climate variables:

fun.Process_climate<-function(df, site, date){
  #make single row dataframes for the processed climate variables:
  # precip
  df.rain<-data.frame(Name = rain_colnames, Value = round(cumsum(df$pr)[rain_cum_and_temp_deltas],2))%>%
    pivot_wider(names_from = Name, values_from = Value)
  # tmax lag
  df.tmax.lag<-data.frame(Name = tmax_lag_colnames, Value = round(df$tmmx[temp_lags+1],2))%>%
    pivot_wider(names_from = Name, values_from = Value)
  # tmax delta
  df.tmax.delta<-data.frame(Name = tmax_delta_colnames, Value = round(df$tmmx[rain_cum_and_temp_deltas+1]-df$tmmx[1],2))%>%
    pivot_wider(names_from = Name, values_from = Value)
  # tmin lag
  df.tmin.lag<-data.frame(Name = tmin_lag_colnames, Value = df$tmmn[temp_lags+1])%>%
    pivot_wider(names_from = Name, values_from = Value)
  # tmin delta
  df.tmin.delta<-data.frame(Name = tmin_delta_colnames, Value = round(df$tmmn[rain_cum_and_temp_deltas+1]-df$tmmn[1],2))%>%
    pivot_wider(names_from = Name, values_from = Value)
  # combine these dataframes and add a site and date column:
  df.date<-bind_cols(df.rain, df.tmax.lag, df.tmax.delta, df.tmin.lag, df.tmin.delta)%>%
    mutate(site_no = vect.NWIS$Name[i], Date = date, .before = 1)
  return(df.date)
}

# normalize from 0 to 1:

normalized<-function(y) {
  
  x<-y[!is.na(y)]
  
  x<-(x - min(x)) / (max(x) - min(x))
  
  y[!is.na(y)]<-x
  
  return(y)
}


add_columns <- function(df, columns){
  new <- rep(NA_character_, length(columns))
  names(new) <- columns
  mutate(df, !!!new)
}



# get watershed percent HSG for a single site:
# *use this function in lapply for multiple sites:
# if function doesnt work (maybe throws error that says file not found), try emptying Users/ryrug/AppData/Local/Temp folder
# 
# site_no<-df.sf.NWIS$Name[1]
# # # 
# sf.df<-df.sf.NWIS

fun.SURRGO_HSG<-function(site_no, sf.df){
  
  print(site_no)
  
  # workflow for one site:
  
  # set function variables:
  
  template<-sf.df$geometry[sf.df$Name==site_no]
  
  # download soil data for site:
  
  l.soils<-FedData::get_ssurgo(template = template, label = site_no) 
  
  # look at map:
  
  mapview(l.soils[[1]])+mapview(template)
  
  # will need to clip but can do that later...
  
  # the mukey in the spatial element can be used as a joining column for the HSG, which is located in l.soils$tabular$muaggatt:
  
  df.sf.soils<-left_join(l.soils[[1]], l.soils$tabular$muaggatt%>%select(mukey, hydgrpdcd)%>%mutate(MUKEY = as.character(mukey), .keep = 'unused'), by = 'MUKEY')
  
  save(df.sf.soils, file = paste0('Processed_Data/SURRGO_dfs/df.sf.soils.', site_no, '.Rdata'))
  
  # load('Processed_Data/df.sf.soils.Rdata')
  
  # lets filter the MU that are NA and replot:
  
  df.sf.soils%>%
    # drop_na(hydgrpdcd)%>%
    mapview(., zcol = 'hydgrpdcd')+mapview(template)
  
  # it looks like the NAs are water!!
  
  # convert to SpatVector to perform spatial analysis:
  
  vect.soils<-vect(df.sf.soils%>%drop_na(hydgrpdcd))
  
  vect.DA.template<-vect(template)
  
  # calcualte watershed area in unit ha:
  
  site_DA.ha<-round(expanse(vect.DA.template, unit="ha", transform=TRUE), 2)
  
  # crop soils vector to draiange area vector:
  
  vect.soils<-terra::crop(vect.soils, vect.DA.template)
  
  # look at map:
  
  plot(x= vect.soils, y='hydgrpdcd')
  
  # add column of each polygons area:
  
  vect.soils$area_ha<-expanse(x= vect.soils, unit = 'ha')
  
  # make a dataframe of the areas and HSG, and calcualte the percent of each HSG:
  
  df.HSG<-data.frame(Area = vect.soils$area_ha, HSG = vect.soils$hydgrpdcd)%>%
    drop_na(HSG)%>%
    group_by(HSG)%>%
    summarize(Watershed_Percent = round(sum(Area)/site_DA.ha, 4))
  
}


#### fun.NHD.RIP.buffers ####

# i<-39

fun.NHD.RIP.buffers<-function(i, df.site_metadata){
  
  print(i)
  
  # set up function inputs:
  
  long<-df.site_metadata$dec_long_va[i]
  lat<-df.site_metadata$dec_lat_va[i]
  
  # workflow for one site:
  
  # run through nhdplusTools workflow:
  
  start_point <- st_sfc(st_point(c(long, lat)), crs = 4269)
  start_comid <- discover_nhdplus_id(start_point)
  
  flowline <- navigate_nldi(list(featureSource = "comid", 
                                 featureID = start_comid), 
                            mode = "upstreamTributaries", 
                            distance_km = 1000)

  subset_file <- tempfile(fileext = ".gpkg")
  subset <- subset_nhdplus(comids = as.integer(flowline$UT$nhdplus_comid),
                           # output_file = subset_file,
                           nhdplus_data = "download",
                           flowline_only = TRUE,
                           # overwrite = TRUE,
                           return_data = TRUE)

  flowline <- subset$NHDFlowline_Network
  # catchment <- subset$CatchmentSP
  # waterbody <- subset$NHDWaterbody
  
  ## Or cando this by reading in from file:
  
  # flowline <- sf::read_sf(subset_file, "NHDFlowline_Network")
  # catchment <- sf::read_sf(subset_file, "CatchmentSP")
  # waterbody <- sf::read_sf(subset_file, "NHDWaterbody")
  
  ## now plot:
  
  # plot(sf::st_geometry(flowline), col = "blue")
  # plot(start_point, cex = 1.5, lwd = 2, col = "red", add = TRUE)
  # plot(sf::st_geometry(catchment), add = TRUE)
  # plot(sf::st_geometry(waterbody), col = rgb(0, 0, 1, alpha = 0.5), add = TRUE)

  # this is so cool!!
  
  # play around with flowline df:
  
  # mapview(flowline, zcol = "ftype")
  
  # my goal with the NHD is to use it to calcualte the landuse in the buffer,
  # so I really only need the StreamRiver ftype. 
  
  flowline1<-flowline%>%filter(ftype == 'StreamRiver')
  
  # create 100 and 800 meter buffer. use st_union to remove overlap:
  
  buffer.100<-st_buffer(flowline1, 100)%>%st_union()
  
  buffer.800<-st_buffer(flowline1, 800)%>%st_union()
  
  # Look at map:
  
  # mapview(buffer.100)+mapview(flowline1)
  
  # looks good
  
  # calculate the watershed 100 and 800 meter percent in major landuse types. to do this:
  
  # Load in the NLCD for the watershed:
  
  rast.NLCD<-rast(names.2006[i])
  
  # reproject to buffer vector data to match raster data:
  
  vect.buffer.100<-vect(buffer.100) 
  vect.buffer.800<-vect(buffer.800)
  
  vect.buffer.100.proj<-terra::project(vect.buffer.100, crs(rast.NLCD))
  vect.buffer.800.proj<-terra::project(vect.buffer.800, crs(rast.NLCD))
  
  # extract frequency tables:
  
  df.buffer.freq.100 <- terra::extract(rast.NLCD, vect.buffer.100.proj, ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2))
  df.buffer.freq.800 <- terra::extract(rast.NLCD, vect.buffer.800.proj, ID=FALSE)%>%group_by_at(1)%>%summarize(Freq=round(n()/nrow(.),2))
  
  # makesure each NLCD class is represented. todo this: left join with legend:
  
  df.buffer.freq.100<-left_join(legend%>%select(Class), df.buffer.freq.100, by = 'Class')%>%replace(is.na(.), 0)
  df.buffer.freq.800<-left_join(legend%>%select(Class), df.buffer.freq.800, by = 'Class')%>%replace(is.na(.), 0)
  
  # reclassify: the GAGES II predictors for NLCD land use are the sum of a few NLCD classes (see the kable table of the variable descriptions made above). To do this:
  
  # create and ordered vector based on legend (legend comes from Ryan_funcitons.R)
  
  Class3.for.G2<-c("WATERNLCD06", "SNOWICENLCD06", rep("DEVNLCD06", 4), "BARRENNLCD06", "DECIDNLCD06", "EVERGRNLCD06", "MIXEDFORNLCD06", NA, "SHRUBNLCD06", "GRASSNLCD06", NA, NA, NA, "PASTURENLCD06", "CROPSNLCD06", "WOODYWETNLCD06", "EMERGWETNLCD06")
  
  # add an identifier to this vector so the column names are slightly different than the GAGES II predictors:
  
  Class3.for.G2<-paste0('R_', Class3.for.G2)
  
  # create a df from this vector (for latter use):
  
  df.Class3<-data.frame(Class = unique(Class3.for.G2)[complete.cases(unique(Class3.for.G2))])
  
  # create new df from the legend dataframe:
  
  legend.for.G2<-legend%>%mutate(Class3 = Class3.for.G2)
  
  # reclassify the NLCD using this new legend and clean up the dataframe:
  
  df.buffer.freq.100<-left_join(df.buffer.freq.100, legend.for.G2%>%select(Class, Class3), by = 'Class')%>%
    mutate(Class = Class3)%>%
    select(-Class3)%>%
    group_by(Class)%>%
    summarise(Freq = sum(Freq))%>%
    pivot_wider(names_from = Class, values_from = Freq)%>%
    mutate(Name = df.sf.NWIS$Name[i], .before = 1)%>%
    mutate(R_FORESTNLCD06 = R_DECIDNLCD06+R_EVERGRNLCD06+R_MIXEDFORNLCD06,
           R_PLANTNLCD06 = R_PASTURENLCD06+R_CROPSNLCD06)%>%
    pivot_longer(cols = -Name)
  
  df.buffer.freq.800<-left_join(df.buffer.freq.800, legend.for.G2%>%select(Class, Class3), by = 'Class')%>%
    mutate(Class = Class3)%>%
    select(-Class3)%>%
    group_by(Class)%>%
    summarise(Freq = sum(Freq))%>%
    pivot_wider(names_from = Class, values_from = Freq)%>%
    mutate(Name = df.sf.NWIS$Name[i], .before = 1)%>%
    mutate(R_FORESTNLCD06 = R_DECIDNLCD06+R_EVERGRNLCD06+R_MIXEDFORNLCD06,
           R_PLANTNLCD06 = R_PASTURENLCD06+R_CROPSNLCD06)%>%
    pivot_longer(cols = -Name)
  
  # I really only need DEV, PLANT, and FOREST.
  # so filtering these, pivoting wider, and merging with G2:
  
  df.buffer.freq.100<-df.buffer.freq.100%>%
    filter(name %in% c('R_DEVNLCD06', 'R_PLANTNLCD06', 'R_FORESTNLCD06'))%>%
    pivot_wider(names_from = name, values_from = value)%>%
    rename_with(~paste0(., "_100"), starts_with("R_"))
  
  df.buffer.freq.800<-df.buffer.freq.800%>%
    filter(name %in% c('R_DEVNLCD06', 'R_PLANTNLCD06', 'R_FORESTNLCD06'))%>%
    pivot_wider(names_from = name, values_from = value)%>%
    rename_with(~paste0(., "_800"), starts_with("R_"))
  
  # create a df of the 100 and 800 meter ripiran predictors:
  
  df.RIP<-left_join(df.buffer.freq.100, df.buffer.freq.800, by = 'Name')
  
  # rename columns to G2, but paste a unique ID to it:
  
  names(df.RIP)<-c('Name', paste0('R_', names(df.G2.reduced%>%select(contains('RIP')))))
  
  return(df.RIP)
}




