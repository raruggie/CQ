# map a folder of csv files as dataframes to the gloabl envrionment:

temp = list.files(pattern="*.csv") 

# then lapply read.csv function to the list of names:

dfs = lapply(temp, read.csv, fileEncoding="UTF-8-BOM")

# add the names of the file names as the names of the dfs in the list, minus the ".csv" (last 4 characters removed)

names(dfs) <- substr(temp,1,nchar(temp)-4) 

# map the items of the list to the global env.:

list2env(dfs,globalenv()) 