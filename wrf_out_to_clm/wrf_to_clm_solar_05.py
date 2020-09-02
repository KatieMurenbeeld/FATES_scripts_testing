# Load modules, insert code, and run programs below #

#module load python/intel/3.5

import xarray as xr
import numpy as np
import pandas as pd
from varlist import var_list
import datetime
import time
import glob

t0 = time.time()

#for year in np.arange(1988, 1989):
wy = 1994
year = 1994
#month = ['01','02','03','04','05','06','07','08','09']
#month = ['10','11','12']
month=['05']

#for y in year:
   #filesy = sorted(glob.glob('home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + 'd02/wrfout_d02_' + str(y) + '-*
for m in month:
    #print(m)
    # Select all files except last (which is a repeat of the first hour of the following water year) for each WY
    #files = sorted(glob.glob('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_d02_' + str(year)+ '-' + str(m) + '-*'))[:-1]
    files = sorted(glob.glob('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(wy) + '/d02/wrfout_d02_' + str(year)+ '-' + str(m) + '-*'))
    print('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_d02_' + str(year) + '-'+ str(m) + '-*')
    print("Total number of files: {}".format(len(files)))

    # open multi-file dataset (this function accepts unix wildcards)
    d = xr.open_mfdataset(files, concat_dim='Time', parallel=True, drop_variables = var_list).chunk({'Time':240})
    
    print(d)
   
    #d = d.load() 
    # Convert times, a character string, to datetimeindex
    
    #time_test = pd.to_datetime(d['Times'][1:]
    #time_test = time_test.strftime("%Y/%m/%d %H:%M:%S")
    #d['Time_test'] = time_test
    #d['Times'] = [datetime.datetime.strptime(x, '%Y/%m/%d %H:%M:%S') for x in d['Times']]
    
    #time = []
    #stimes = vars(d['Times'])
    #print(d['Times'])
    #for t in stimes:
    #    time.append(pd.to_datetime(t, errors='raise', format='"%Y-%m-%d_%H:%M:%S"'))
    #    print(time)
    
    index = 48
    print(index)
    print(d['XTIME'])
    
    d = d.load()

    d['XTIME'][index] = d['XTIME'][index -1] + np.timedelta64(1, 'h')

    #d['Times'].astype('str') 
    #d['Times'] = ('Time', pd.to_datetime(d['Times'].values, format='"%y-%m-%d_%H:%M:%S"'))
    
    # Swap time and Times
    d = d.swap_dims({'Time':'XTIME'})
    #d = d.swap_dims({'Time':'Times'})    

    # Set coordinates for the new xarray dataset
    d3_xlat = np.array(d['XLAT']) # Create a 3-D array from the XLAT
    d2_xlat = d3_xlat[0,:,:] # Make new 2-D array for XLAT
    d3_xlon = np.array(d['XLONG']) # Create a 3-D array from the XLON
    d2_xlon = d3_xlon[0,:,:] # Make a new 2-D array for XLON


    # Get 1-hourly data for the desired variables
    new_array = d[['SWNORM']].resample(XTIME = '1H').mean(dim = 'XTIME') # create daily means of few variables
    #print('XTIME')
    #new_array.rename({'SWNORM':'FSDS'}) # rename to be consistent with CLM(FATES)

    # Adjust the meta data
    #new_array['SWNORM'].attrs = [('description','1-HOURLY MEAN NORMAL SHORT WAVE FLUX AT GROUND SURFACE (SLOPE-DEPENDENT)'), ('units','W m^2')]

    # Create a new xarray dataset
    ds = xr.Dataset({'FSDS': (['time','x','y'], new_array['SWNORM']),
                'LATIXY': (['x','y'], d2_xlat),
                'LONGXY': (['x','y'], d2_xlon)},
                coords = {'lat': (['x','y'], d2_xlat),
                        'lon': (['x','y'], d2_xlon),
                        'time': np.array(d['XTIME'])})
    ds['FSDS'].attrs = [('description','1-HOURLY MEAN NORMAL SHORT WAVE FLUX AT GROUND SURFACE (SLOPE-DEPENDENT)'), ('units','W m^2')]

    # Write new netcdf file
    print("Writing output to disk")
    ds.to_netcdf("/home/kmurenbeeld/scratch/wrf_to_clm/clmforc_WRF30d02_c2018_1kmx1km_" + str(year) + "_" + str(m) + "_Solar.nc")

t1 = time.time()
print("Total time to create this subset was:", t1 - t0, "seconds.")

