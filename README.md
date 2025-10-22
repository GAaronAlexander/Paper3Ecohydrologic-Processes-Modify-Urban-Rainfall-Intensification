This repo contains the WRF model used in the study *Ecohydrologic Processes Modify Urban Rainfall Intensification in Land-Atmosphere Simulations*. The folders are structured as follows: 

## WRF 
This folder contains the mother directory of the WRF model *(WRF-wrf-v-4.4.1-noahhue)* and a directory with example namelists *(namelist-example)* to run WRF in this study. 

To complie WRF, one needs to follow standard instructions found in the [WRF Tutorial] (https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php) with one excepetion: 
- When configuring the model run: 
	```
	./configure hue
	```
	This will 'hide' some variables to allow for mosaic and HUE variables to be active. Note that compiling Noah-MP HUE will remove the ability to run the CLM model with the same executable. 

Within the WRF folder is also the MPTABLE.MPTABLE.HUE.7AUG2025.TBL that contains new parameters for this study. Note that we only employ landcovers 42 and 43, as well as 41 and 47 for this study. 

## Observed-data
This folder contains the observed data used in *Ecohydrologic Processes Modify Urban Rainfall Intensification in Land-Atmosphere Simulations*.

### AMDAR
Aircraft Meteorological Data Relay (AMDAR) is hosted by the World Meteorological Organization. This folder contains: 
1. Download-script: contains a python script to download AMDAR data in its raw from 
2. CleaningAMDAR: A script that outlines how data was pooled and averaged. This script does require surface meteorological observations, which were provided in folder KMKE
3. The resulting files, which contain temprature, humidity, U, and V winds, in netcdf format, with standard deviations!

### MMSD-rainfall
This study does use raingage data that was provided by the Milwaukee Metropolitan Sewerage District. The full dataset is not public, but the data used for this study is available for comparisons. Note that this data is organized in netcdf files, with 20 gages and 5 minute time periods of measurement. 

The 20 gages are **not** QAQC'd in these files, meaning that for any given storm, one or more stations could be measuring nan's (e.g. offline). The stations that were inactive does change between each storm, so one should remove nans carefully!

