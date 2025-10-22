import xarray as xr
import numpy as np
import pandas as pd
from scipy import interpolate
from statsmodels.nonparametric.smoothers_lowess import lowess as  sm_lowess
#from pykrige import OrdinaryKriging
#import gstools as gs
import os 
import matplotlib.pyplot as plt
import warnings
from datetime import datetime


warnings.filterwarnings('ignore') #there are crazy warnings that are going to be printed out due to nans
from scipy.stats.mstats import winsorize

def extract_pandas_data(dt_s,dt_e,pandas_data_location):
    '''
    reads in data from mesowest, grabs the columsn that are needed
    then resamples to an hourly dataset
    finally,subsets the data

    ensure that you change the use-cols list. I usually change the headers to be easier
    '''
    dataset = pd.read_csv(pandas_data_location,index_col=['time'],
        usecols=['time','air_temp','wind_speed','wind_direction','pressure'],
        header=0,parse_dates=True)
    dataset_hourly = dataset.resample('H').first()
    dataset_hourly_sub = dataset_hourly[dt_s:dt_e]
    dataset_hourly_sub = dataset_hourly_sub.tz_localize(None)
    return(dataset_hourly_sub)

def find_profile_numbers(air_data,code):
    '''
    finds the profile locations in the netcdf data
    '''
    locs = np.argwhere(air_data['profileAirport'].values.astype('str') == code)

    return(locs)
def change_altitude_to_pressure(air_data,locs,j):
    '''
    takes the altimiter pressure 
    and changes to a pressure that is referenced by 1013.25 mb
    note: air_data altitude is in meters, but if not, 
    change Kf to .3048
    '''
    ## following the WMO pressure data

    Kf = 1

    interior = (1 - (6.5*10**(-3))*(Kf*air_data.altitude[locs[j]].values)/288.5)

    pressure = 1013.25 * (interior**5.2559)
    return(pressure)
def hypsometric_eq(Tv,P1,P2):
    '''
    Takes the hypsometric equation to find the delta h between locations
    Note: we are assuming that T = Tv (e.g. there is very little water 
    in the air)

    T is temperature in K
    P1 & P2 are pressures in hpa
    '''
    R = 287 #J/K kg
    gravity = 9.81 #m/s^2

    # h = z2 = z1
    h = (R*Tv/gravity) * np.log(P1/P2)
    return(h)
def barometric_eq(Tb, Pb,Lb,hb,T,h):
    '''
    Takes heights and temperatures from airplanes, and 
    the: 
    Tb weather station temperature K
    Pb is weather pressure hpa
    Lb is the local lapse rate in K/m
    hb is the height of the weather station above sea level
    T is temperature from airplanes K
    h is heights in meters
    '''
    # calculate regression
    cons = (0.02896*9.81)/(8.3143*-Lb)
    # calcualte pressure
    P = Pb*((Tb - ((h+hb))*-Lb)/Tb)**cons
    # return that pressure
    return(P)
def get_U_V(ws,wd):
    '''
    given a wind speed and wind direction measurement
    derive the U and V vectors
    ws is meters per second
    wd is wind in degrees
    '''
    U = ws*np.cos(((270-wd)%360)*np.pi/180)
    V = ws*np.sin(((270-wd)%360)*np.pi/180)


    return(U,V)
def create_xarray_dataset():
    '''
    takes the average regressed data 
    and creates a netcdf!
    '''
    return(2)

def haversine(coord1: object, coord2: object):
    
  
    # Coordinates in decimal degrees (e.g. 2.89078, 12.79797)
    lat1,lon1 = coord1
    lat2,lon2 = coord2
    

    R = 6371000  # radius of Earth in meters
    phi_1 = np.radians(lat1)
    phi_2 = np.radians(lat2)

    delta_phi = np.radians(lat2 - lat1)
    delta_lambda = np.radians(lon2 - lon1)

    a = np.sin(delta_phi / 2.0) ** 2 + np.cos(phi_1) * np.cos(phi_2) * np.sin(delta_lambda / 2.0) ** 2
    
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    meters = R * c  # output distance in meters
    km = meters / 1000.0  # output distance in kilometers

    return(km)

def kriging_fit(heights,variable,heights_wanted,nuggets,length_scale):

                
                #first fit the optimal amount of bins and the gamma
            
                bin_center, gamma = gs.vario_estimate((heights[~np.isnan(variable)], np.zeros(heights[~np.isnan(variable)].shape)), variable[~np.isnan(variable)])

                #then we are going to fit a variogram model
                fit_model = gs.Cubic(variance = 12,nugget=nuggets,len_scale=length_scale,dim=1)
                para, pcov, r2, = fit_model.fit_variogram(bin_center,  gamma, return_r2=True,len_scale=False)
                
                
                ok = gs.krige.Ordinary(
                    model=fit_model,
                    cond_pos=heights,
                    cond_val=variable,
                    
                    
)
               
             
                ok((heights_wanted))
                y_pred = ok.field
                y_std = np.sqrt(ok.krige_var)
       
                return(y_pred,y_std)

def lowess_and_confidence_intervals(heights_,variable_,heights_want):
    
    heights_inorder, variable_inorder = sm_lowess(variable_, heights_,  frac=1./10., 
                           it=5, return_sorted = True).T
   
  
    variable_commongrid = interpolate.interp1d(heights_inorder, variable_inorder, 
                                        fill_value='extrapolate')(heights_want)
    
    
    ## now get confidence intervals
    
    samples = np.random.choice(len(heights_inorder), len(heights_inorder), replace=True)
    
    K = 1000
    smooths = np.zeros((heights_want.shape[0],K))
    for g in range(K):
        smooths[:,g] =smooth(heights_, variable_, heights_want)
    

    
    return(variable_commongrid,smooths)


def smooth(heights, variable, heights_want):
    
    samples = np.random.choice(len(heights), len(heights), replace=True)
    variable_s = variable[samples]
    heights_s = heights[samples]
    y_sm = sm_lowess(variable_s,heights_s, frac=1./10., it=2,
                     return_sorted = False)
    # regularly sample it onto the grid
    random_sample_smooth = interpolate.interp1d(heights_s, y_sm, 
                                        fill_value='extrapolate')(heights_want)
    return random_sample_smooth

def clean_by_heights(heights_original, variable_origional, HEIGHTS_=np.array([250,1000])):
    
    # this function takes in a height variable and a co-located variable (Temperature, winds, etc) and applies a Median Absolute 
    # deviation threshold of 90%. This is done over the HEIGHTS_ variable
    
    HEIGHTS_loop = np.append(HEIGHTS_, HEIGHTS_[-1])
    
    heights_final = []
    variable_final = []
    
    for i in range(len(HEIGHTS_loop)):

    
        if i == 0:
            bool_loop_heights = np.where(heights_original <= HEIGHTS_loop[i],True,False)
        elif i == 2:
            bool_loop_heights = np.where(heights_original > HEIGHTS_loop[i],True,False)
            
        else:
            bool_loop_heights = np.where((heights_original > HEIGHTS_loop[i-1]) & (heights_original <= HEIGHTS_loop[i]),True,False)
        
        
        x = np.abs((variable_origional[bool_loop_heights]) - np.nanmedian(variable_origional[bool_loop_heights]))
        quants = np.nanpercentile(x, 90)
        heights_bool_MAD = np.where(x<quants,True,False)
        
        heights_final.append(heights_original[bool_loop_heights][heights_bool_MAD])
        variable_final.append(variable_origional[bool_loop_heights][heights_bool_MAD])
    
    
    return (np.concatenate(heights_final),np.concatenate(variable_final))


# we can ignore
#%%

start_dates = ['2019-09-10 19:00','2016-09-06 18:00','2019-07-19 20:00','2020-07-08 18:00','2021-08-06 01:00','2018-06-17 19:00',
               '2015-07-06 00:00','2018-09-02 13:00','2018-09-04 16:00','2017-07-11 08:00','2014-06-17 00:00']

end_dates =   ['2019-09-11 19:00','2016-09-07 18:00','2019-07-20 20:00','2020-07-09 18:00','2021-08-07 01:00','2018-06-18 19:00',
               '2015-07-07 00:00','2018-09-03 13:00','2018-09-05 16:00','2017-07-12 08:00','2014-06-18 00:00']

wrf_STARTS =  ['2019-09-08 12:00','2016-09-05 00:00','2019-07-18 00:00','2020-07-06 00:00','2021-08-05 12:00','2018-06-16 00:00',
               '2015-07-04 00:00','2018-09-01 00:00','2018-09-01 00:00','2017-07-09 00:00','2014-06-15 00:00']
wrf_ENDS =    ['2019-09-13 00:00','2016-09-08 12:00','2019-07-21 12:00','2020-07-11 00:00','2021-08-09 00:00','2018-06-21 00:00',
               '2015-07-08 00:00','2018-09-06 12:00','2018-09-06 12:00','2017-07-13 00:00','2014-06-19 00:00']

wrf_save_name = ['12SEP2019','07SEP2016','20JULY2019','09JULY2020','08AUG2021','18JUNE2018','06JULY2015','03SEP2018','05SEP2018','12JULY2017','18JUNE2014']


location_of_interest = f'MKE'
final_data_name_prefix = f'{location_of_interest}_AMDAR_LOWESS_NONSEABREEZE_WS_removeoutliers.nc' 
save_location = f'./'
pandas_data_location = f'../../KMKE/'
data_location_raw = f'/Volumes/Backup-Plus/dissertation_chapter_3/AMDAR_data/' # put this where you have saved the outputs
height_wanted = 2000 # in meters, the number up to which you want data output!
height_of_airport = 188 #m




## get list of days
nonseabreeze = [0,1,2,4,5,6,7,8,10] #these are the non-seabreeze events
seabreeze = [3,9]

files_to_load = []
pandas_data = []

for event in nonseabreeze:
    starts_ = start_dates[event]
    ends_ = end_dates[event]
    date_range = pd.date_range(start=starts_,end=ends_,freq='H')
    
    
    date_range_bool = np.where(((date_range.hour<23) & (date_range.hour>=17) & (date_range.day != ends_[8:10])),True,False)
    date_range = date_range[date_range_bool]
    
    
    single_events = extract_pandas_data(starts_,ends_,pandas_data_location+f'KMKE_{wrf_save_name[event]}.csv')
    print(single_events.shape)
    print(single_events)
    if single_events.shape[0] == date_range_bool.shape[0]:
        pandas_data.append(single_events[date_range_bool])
    elif single_events.shape[0] > date_range_bool.shape[0]:
        pandas_data.append(single_events[date_range_bool[:single_events.shape[0]-1]])
    else:
        date_range_bool = date_range_bool[0:-1]
        pandas_data.append(single_events[date_range_bool[:single_events.shape[0]]])
    
    
    
    for date__ in date_range:
        files_to_load.append(f'{data_location_raw}{wrf_save_name[event]}/{str(date__.year)}{str(date__.month).zfill(2)}{str(date__.day).zfill(2)}_{str(date__.hour).zfill(2)}00')
        
        


pandas_all = pd.concat(pandas_data)
surface_data = pandas_all

    
    








## set up daily variables and things that do not need to be renewed in a for loop 

# setting up the array of heights, incremented in 20 meters
new_heights = np.arange(0,height_wanted,50)
# the height of 2 meters is a constant and is always defined
new_heights[0] = 2
# how big should the means be?
height_shape = new_heights.shape


temp_all = []
heights_all = []
U_all = []
V_all = []
Tdew_all = []

n = 0

for d,files in enumerate(files_to_load):

    v = xr.open_dataset(files,engine='netcdf4')




    

    locations = find_profile_numbers(v,location_of_interest)

    print(locations)
    
    n+=locations.shape[0]
    print(n)

    # get the number of profiles
    shape_of_temp = locations.shape[0]

    if shape_of_temp == 0:
        continue

    else:
        #if there are actually things to be able to averages


        pres_obs = surface_data.iloc[d].pressure/100
        temp_obs = surface_data.iloc[d].air_temp+273.15
        
        ws_obs = surface_data.iloc[d].wind_speed
        wd_obs = surface_data.iloc[d].wind_direction

       

        h2 = np.zeros(200) # heights


        latitude_ = []
        longitude_ = []


        for j in range(0,shape_of_temp):
            
            pres = change_altitude_to_pressure(v,locations,j)

            Tavg = (temp_obs + v.temperature[locations[j]][0,0].values)/2

            # create hypsometric equation for first layer
            h = hypsometric_eq(Tavg,pres[0,0],pres_obs) 


            if np.abs(h) > 200: #if above or below 200 meters, discard
                
                continue

            else: # fill in the other heights
                h2[0] = h 
                
                t_avg2 = (v.temperature[locations[j]][0,1:].values + v.temperature[locations[j]][0,:-1].values)/2

                h_ = hypsometric_eq(t_avg2,pres[0,:-1],pres[0,1:]) #
                h2[1:] = h_


            heights = np.cumsum(h2) #cumulative sum the heights


            # derive the local lapse rate for the barometric pressure equation
            a = v.temperature[locations[j]][0,1:].values-v.temperature[locations[j]][0,:-1].values
            b = heights[1:]-heights[:-1]


            # calcualte the pressure for temperature
            press= barometric_eq(temp_obs,pres_obs,a/b,height_of_airport,v.temperature[locations[j]][0,:-1].values,heights[:-1])
            pot_temp = v.temperature[locations[j]][0,:-1].values* (1013/press)**.286
            
            U_raw,V_raw = get_U_V(v.windSpeed[locations[j]][0,:-1].values,v.windDir[locations[j]][0,:-1].values)
            
            qvapor_sat = 6.112 * np.exp (17.62 * (v.temperature[locations[j]][0,:-1].values.flatten() - 273.15)/(243.12 + (v.temperature[locations[j]][0,:-1].values.flatten() - 273.15)))
            qvapor_temp = qvapor_sat * v.relHumidity[locations[j]][0,:-1].values.flatten()/100
            qvapor_temp = 0.622*qvapor_temp/(press - qvapor_temp)
            qvapor_final = qvapor_temp/(qvapor_temp+1)
            
            

            #filter heights and all other values 
            heights_ = heights[:-1].flatten()
            temp_all_ = pot_temp.flatten()
            U_all_ = U_raw.flatten()
            V_all_ = V_raw.flatten()
            Tdew_all_ = qvapor_final
            # remove nans based on heights_
            filters = ~np.isnan(heights_)
            
            #nab the lat and lon
            lat_ = v.trackLat[locations[j]][0,1:].values
            lon_ = v.trackLon[locations[j]][0,1:].values
            
           
            
            distance = haversine([lat_,lon_], [np.ones(lat_.shape)*(v.latitude[locations[j]].values[0]),np.ones(lat_.shape)*(v.longitude[locations[j]].values[0] )])
            #x = haversine([51.510357,-0.116773, ], [ 38.889931,-77.009003])

            distance_bool = np.where((distance[filters]>=15) & (heights_[filters]<=1500) ,False,True)
           
            # we are trying to pool all our data!
            heights_all.append(heights_[filters][distance_bool])
            temp_all.append(temp_all_[filters][distance_bool])
            U_all.append(U_all_[filters][distance_bool])
            V_all.append(V_all_[filters][distance_bool])
            Tdew_all.append(Tdew_all_[filters][distance_bool])


# flatten the arrays

heights_all2=  np.concatenate(heights_all).ravel()
temp_all2=  np.concatenate(temp_all).ravel()
U_all2=  np.concatenate(U_all).ravel()
V_all2=  np.concatenate(V_all).ravel()
Tdew_all2=  np.concatenate(Tdew_all).ravel()

# sort the arrays
idx = np.argsort(heights_all2)

heights_all2 = np.take_along_axis(heights_all2,idx,axis=0)
temp_all2 = np.take_along_axis(temp_all2.flatten(),idx,axis=0)
Tdew_all2 = np.take_along_axis(Tdew_all2,idx,axis=0)
U_all2 = np.take_along_axis(U_all2,idx,axis=0)
V_all2 = np.take_along_axis(V_all2,idx,axis=0)


# filter heights to speed up regression
filter_heights = np.where((heights_all2 <= height_wanted) & (heights_all2 >= -20),True,False)

heights_below = heights_all2[filter_heights]
heights_below = np.where( heights_below < 0,np.abs(heights_below),heights_below)
    


temp_below = temp_all2[filter_heights]
Tdew_below = Tdew_all2[filter_heights]
U_below = U_all2[filter_heights]
V_below = V_all2[filter_heights]


        





print('regressing and bootstrapping - Temperature')
heights_cleaned, theta_cleaned = clean_by_heights(heights_below, temp_below)

temperture_regressed,temperature_samples = lowess_and_confidence_intervals(heights_cleaned,theta_cleaned,new_heights)

temperature_05 =  np.nanquantile(temperature_samples,0.05,axis=1) 

temperature_95 =  np.nanquantile(temperature_samples,0.95,axis=1) 

print('regressing and bootstrapping - U-Wind')
heights_cleaned, U_cleaned = clean_by_heights(heights_below, U_below)

U_regressed,U_samples = lowess_and_confidence_intervals(heights_cleaned,U_cleaned,new_heights)

U_05 =  np.nanquantile(U_samples,0.05,axis=1) 

U_95 =  np.nanquantile(U_samples,0.95,axis=1) 

print('regressing and bootstrapping - V-Wind')
heights_cleaned, V_cleaned = clean_by_heights(heights_below, V_below)

V_regressed,V_samples = lowess_and_confidence_intervals(heights_cleaned,V_cleaned,new_heights)

V_05 =  np.nanquantile(V_samples,0.05,axis=1) 

V_95 =  np.nanquantile(V_samples,0.95,axis=1) 

print('regressing and bootstrapping - Water Vapor')
heights_cleaned, Q_cleaned = clean_by_heights(heights_below, Tdew_below)

Q_regressed,Q_samples = lowess_and_confidence_intervals(heights_cleaned,Q_cleaned,new_heights)

Q_05 =  np.nanquantile(Q_samples,0.05,axis=1) 

Q_95 =  np.nanquantile(Q_samples,0.95,axis=1) 



dataset_name = f'{final_data_name_prefix}'


pot_temp_mean_final = xr.DataArray(
                data = temperture_regressed,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='Mean Potential Temperature Profile Above Ground Level',
                        units='K',
                        standard_name = 'air_potential_temperature',
                        long_name='Mean Potential Temperature',
                        missing_value=np.nan
                ),
)


pot_temp_05 = xr.DataArray(
                data = temperature_05,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='5th Quantile Potential Temperature Profile Above Ground Level',
                        units='K',
                        standard_name = 'air_potential_temperature_q05',
                        long_name='Q05 Potential Temperature',
                        missing_value=np.nan
                ),
)

pot_temp_95 = xr.DataArray(
                data = temperature_95,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='95th Quantile Potential Temperature Profile Above Ground Level',
                        units='K',
                        standard_name = 'air_potential_temperature_q95',
                        long_name='Q95 Potential Temperature',
                        missing_value=np.nan
                ),
)


q_mean_final = xr.DataArray(
                data = Q_regressed,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='Mean Specific Humidity Profile Above Ground Level',
                        units='kg kg',
                        standard_name = 'specific_humidity',
                        long_name='Mean Specific Humidity',
                        missing_value=np.nan
                ),
)


q_05 = xr.DataArray(
                data = Q_05,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='5th Quantile Specific Humidity Profile Above Ground Level',
                        units='kg kg',
                        standard_name = 'specific_humidity_q05',
                        long_name='Q05 Specific Humidity',
                        missing_value=np.nan
                ),
)

q_95 = xr.DataArray(
                data = Q_95,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='95th Quantile Specific Humidity Profile Above Ground Level',
                        units='kg kg ',
                        standard_name = 'specific_humidity_q95',
                        long_name='Q95 pecific Humidity',
                        missing_value=np.nan
                ),
)


U_mean_final = xr.DataArray(
                data = U_regressed,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='Mean U Wind Profile Above Ground Level',
                        units='m s-1',
                        standard_name = 'U_wind',
                        long_name='Mean U Wind',
                        missing_value=np.nan
                ),
)

U_05 = xr.DataArray(
                data = U_05,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                       description='5th Quantile U Wind Profile Above Ground Level',
                       units='m s-1',
                       standard_name = 'U_wind_q05',
                       long_name='Q05 U Wind',
                       missing_value=np.nan
                ),
)

U_95 = xr.DataArray(
                data = U_95,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='95th Quantile U Wind Profile Above Ground Level',
                        units='m s-1',
                        standard_name = 'U_wind_q95',
                        long_name='Q95 U Wind',
                        missing_value=np.nan
                ),
)



V_mean_final = xr.DataArray(
                data = V_regressed,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='Mean V Wind Profile Above Ground Level',
                        units='m s-1',
                        standard_name = 'V_wind',
                        long_name='Mean V Wind',
                        missing_value=np.nan
                ),
)


V_05 = xr.DataArray(
                data = V_05,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                       description='5th Quantile V Wind Profile Above Ground Level',
                       units='m s-1',
                       standard_name = 'V_wind_q05',
                       long_name='Q05 V Wind',
                       missing_value=np.nan
                ),
)

V_95 = xr.DataArray(
                data = V_95,
                dims=["height_abl"],
                coords=dict(
                        height_abl=(["height_abl"],new_heights),
                        
                ),
                attrs=dict(
                        description='95th Quantile V Wind Profile Above Ground Level',
                        units='m s-1',
                        standard_name = 'V_wind_q95',
                        long_name='Q95 V Wind',
                        missing_value=np.nan
                ),
)





dataset_final = xr.Dataset(
            data_vars=dict(
                Pot_Temp_mean=pot_temp_mean_final,
                
                Qvapor_mean=q_mean_final,
                    
                U_mean=U_mean_final,
                    
                V_mean=V_mean_final,
                
                Pot_Temp_q05=pot_temp_05,
                
                Qvapor_q05=Q_05,
                    
                U_Q05=U_05,
                    
                V_Q05=V_05,
                
                Pot_Temp_q95=pot_temp_95,
                
                Qvapor_q95=Q_95,
                    
                U_Q95=U_95,
                    
                V_Q95=V_95,
                

                ),
            coords=dict(
                height_abl=(["height_abl"],new_heights),
                ),
            attrs=dict(
                title='AMDAR Profiles',
                summary=f'AMDAR profiles for {location_of_interest} Airport, generated by pooling data and LOWESS to gain  \'mean\' profiles and quantify uncertainity via bootstrapping 1000 times',
                date_created=f'{datetime.now()}',
                creator_name='G. Aaron Alexander',
                creator_email='gaalexander3@wisc.edu'
                ),


    )

dataset_final.height_abl.attrs['standard_name'] = 'height_abl'  # Optional
dataset_final.height_abl.attrs['long_name'] = 'height above ground level'
dataset_final.height_abl.attrs['units'] = 'm'

print(f'writing {dataset_name}')
dataset_final.to_netcdf(f'{save_location}{dataset_name}')

    


# q = hypsometric_eq(315,1000,pres)

