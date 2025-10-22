import urllib.request
import time
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
import pandas as pd
import os
'''
script created by G. Aaron Alexander on 10 may 2023
intent to download AMDAR data in parallel 
via urllib.request and multiprocessing!

note that this is embarrisingy parallel, 
so there is no need to think about appending data 

'''


directory = '/Volumes/Backup-Plus/dissertation_chapter_3/AMDAR_data/20JULY2019/' # where are we dowloading data
dates = pd.date_range(start='2019-07-18 00:00',end='2019-07-21 00:00',freq='h') # date range of data to download


# make the direcotry
os.system(f'mkdir {directory}')

#create iterables
urls = []
filenames = []

def download_url(args):
	## times, unzips, requests via a try except loop
    t0 = time.time()
    url, fn = args[0], args[1]
    try:
    	urllib.request.urlretrieve(url,fn)
    	# note that amdar data is zipped, so we need to unzip 
    	os.system(f'gunzip -f {fn}')
    	return(url,time.time() - t0)
    except: 
        print(f'URL not found: {url}')

## mapping to a parallel process, and then run it
def download_parallel(args):
    cpus = cpu_count()
    results = ThreadPool(cpus - 1).imap_unordered(download_url, args)
    for result in results:
    	# printing instead of a 
    	if result is not None:
       		print('url:', result[0], 'time (s):', result[1])

for i,d in enumerate(dates):
	# create the file name
	filename_1 = f'{d.year}{str(d.month).zfill(2)}{str(d.day).zfill(2)}_{str(d.hour).zfill(2)}00.gz'
	# create the location you are downlowading to 
	filename_2 = f'{directory}/{d.year}{str(d.month).zfill(2)}{str(d.day).zfill(2)}_{str(d.hour).zfill(2)}00.gz'
	# create the URL location
	url_ = f'https://madis-data.cprk.ncep.noaa.gov/madisPublic1/data/archive/{str(d.year)}/{str(d.month).zfill(2)}/{str(d.day).zfill(2)}/point/acarsProfiles/netcdf/{filename_1}'
	# append this data
	urls.append(url_)
	filenames.append(filename_2)

# zip the data
inputs = zip(urls,filenames)

#download the data
download_parallel(inputs)
