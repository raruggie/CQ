Log of what is done in each code file:

NWIS_XXX.R files:

The NWIS_XXX.R files are replicates of each other but the consitient is swapped out. 
They also varry slightly if something was found when going through the particuale consituent

NWIS.R (XXX refers to any of the nutrients or just the plain NWIS.R which is TP):
- NWIS Query
- Building CQ df (merging flow and nutrient dataframes from query)
- tradeoff matrix (number of samples as a funciton of watershed size)
- Apply filters to number of sites (20 paired CQ observations, after 2001, long island, GAGES 2)
- Fitting breakpoints to CQ relaitonship
- calcualting average annual consituent Yield
- Correlaiton analysis
- Categorzing land use (adjusted from USGS thresholds)
- grouping CQ curves (staitonary, mobilizing, diltuionary, complex)
- Triangle plot - mobilizing strength dependent on disturbance (Ag and devloped on axes)
- MLR
The addiitonal NWIS.R files, which have the consituent name afer (e.g. NWIS_N03.R) have these as well

NWIS_ALL.R file: this is used to take a look at the 17 sites that overlapped between TP, TN, and SRP

NWIS_Datalayers_Download.R file:
- Download an early and recent CDL for a df of watershed shapefiles
- regular and reclassified analysis of the CDL
- heatmap of the difference between early and recent CDLs
- recreating the GAGES II predictors: 
 - land use: using plain CDL for GAGES II CDL and reclassified CDL for GAGES NLCD 06
 - elevation: NED

NWIS_Delineate_and_WWTP.R file:
- delinate watersheds using SS API
- calculated delineation error and determing which sites did not delineate correctly 
- import shapefiles of sites that did not work with API and had todo by hand on SS website
- finding which WWTP (point shapefile) in watershed (polygon shapefiles)

FingerLakesPresentation.R file:
- Map of sites, CQ curves, Triangle plots, are all duplicates of the code in NWIS.R (and the other consituents)
- Combining MLR (turing the MLR from each consituent code into a single table of multiple consituents AND CQ parameters (intercept and AANY))
- kable tables of site result list (in manuscript) and MLR model comparison
- Conceptual CQ diagram
- test(s) if elevaiton is a proxy for land use

