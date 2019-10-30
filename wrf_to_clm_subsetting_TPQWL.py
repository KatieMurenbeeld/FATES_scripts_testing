# Load modules, insert code, and run programs below #

#module load python/intel/3.5

import xarray as xr
import numpy as np
from varlist import var_list
import time
import glob

t0 = time.time()


#for year in np.arange(1988, 1989):
year = 1988
#month = ['01','02','03','04','05','06','07','08','09','10','11','12']
month = ['03','04']
#for y in year:
   #filesy = sorted(glob.glob('home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + 'd02/wrfout_d02_' + str(y) + '-*
for m in month:
    #print(m)
    # Select all files except last (which is a repeat of the first hour of the following water year) for each WY

    #files = sorted(glob.glob('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_d02_' + str(year)+ '-' + str(m) + '-*'))[:-1]
    files = sorted(glob.glob('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_d02_' + str(year)+ '-' + str(m) + '-*'))
    #files = sorted(glob.glob('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_*'))[:-1]
    print('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_d02_' + str(year) + '-'+ str(m) + '-*')
    print("Total number of files: {}".format(len(files)))

    # open multi-file dataset (this function accepts unix wildcards)
    d = xr.open_mfdataset(files, concat_dim='Time', parallel=True).chunk({'Time':240})
    print(d)
    # Swap time and XTIME
    d = d.swap_dims({'Time':'XTIME'})

    # Calculate wind speed from the U and V components
    d['WIND'] = np.sqrt(d['U10']**2 + d['V10']**2)

    # Set coordinates for the new xarray dataset
    d3_xlat = np.array(d['XLAT']) # Create a 3-D array from the XLAT
    d2_xlat = d3_xlat[1,:,:] # Make new 2-D array for XLAT
    d3_xlon = np.array(d['XLONG']) # Create a 3-D array from the XLON
    d2_xlon = d3_xlon[1,:,:] # Make a new 2-D array for XLON

    # Get 1-hourly data for the desired variables
    new_array = d[['T2','Q2','PSFC','SWNORM','GLW','WIND']].resample(XTIME = '1H').mean(dim = 'XTIME') # create daily means of few variables
	new_array = new_array.rename({'T2':'TBOT','Q2':'QBOT','PSFC':'PSRF','GLW':'FLDS'}) # rename to be consistent with CLM(FATES)

    # Create a new dataset
    ds = xr.Dataset({'TBOT':(['time','x','y'], new_array['TBOT']),
                      'QBOT': (['time','x','y'], new_array['QBOT']),
                      'PSRF': (['time','x','y'], new_array['PSRF']),
                      'FLDS': (['time','x','y'], new_array['FLDS']),
                      'WIND': (['time','x','y'], new_array['WIND'])},
                {'LATIXY': (['x','y'], d2_xlat),
                'LONGXY': (['x','y'], d2_xlon)},
                coords = {'lat': (['x','y'], d2_xlat),
                        'lon': (['x','y'], d2_xlon),
                        'time': np.array(d['XTIME'])})
    # Adjust the meta data
	ds['TBOT'].attrs = [('description','1-HOURLY MEAN GRID SCALE TEMPERATUTE'), ('units','K')]
	ds['QBOT'].attrs = [('description','1-HOURLY MEAN GRID SCALE SPECIFIC HUMIDITY'), ('units','kg/kg')]
	ds['PSRF'].attrs = [('description','1-HOURLY MEAN SURFACE PRESSURE'), ('units','Pa')]
	ds['FLDS'].attrs = [('description','1-HOURLY MEAN LONG WAVE FLUX AT GROUNF SURFACE'),('units','Wm^2')]
	ds['WIND'].attrs = [('description','1-HOURLY MAXIMUM GRID SCALE TEMPERATURE'), ('units','K')]

    # Write new netcdf file
    print("Writing output to disk")
    ds[['TBOT','QBOT','PSRF','WIND','FLDS']].to_netcdf("/home/kmurenbeeld/scratch/wrf_to_clm/clmforc_WRF30d02_c2018_1kmx1km_" + str(year) + "_" + str(m) +  "_TPQWL.nc")

      t1 = time.time()
print("Total time to create this subset was:", t1 - t0, "seconds.")
