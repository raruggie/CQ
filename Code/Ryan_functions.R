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

# parameter_code<-'00665' # pcode for TP
# 
# df.flow_query<-df.NWIS.Q_sites

fun.df.Pair_consit_flow<-function(parameter_code, df.flow_query, n_samples = 20) {
  
  # download metadata for the sites in NY with data of the constituent of interest:
  
  df.consit_fun<-whatNWISdata(stateCd="NY",parameterCd=parameter_code)
  
  # reduce duplicate site entries:
  
  df.consit_fun<-df.consit_fun%>% 
    group_by(site_no)%>% 
    slice(which.max(count_nu)) 
  
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


