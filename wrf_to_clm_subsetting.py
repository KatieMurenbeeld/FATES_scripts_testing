# Load modules, insert code, and run programs below #

#module load python/intel/3.5

import xarray as xr
import numpy as np
from varlist import var_list
import time
import glob

t0 = time.time()


# Decummulate Precipitation
def calc_precip(cum_precip, bucket_precip):

    total_precip = cum_precip + bucket_precip * 100.0
    PRCP = np.zeros(total_precip.shape)

    for i in np.arange(1,PRCP.shape[0]):

        PRCP[i,:,:] = total_precip[i,:,:].data - total_precip[i-1,:,:].data

    return PRCP

#for year in np.arange(1988, 1989):
year = 1988
#month = ['01','02','03','04','05','06','07','08','09','10','11','12']
month = ['03','04']
#for y in year:
   #filesy = sorted(glob.glob('home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + 'd02/wrfout_d02_' + str(y) + '-*
for m in month:
    #print(m)
    # Select all files except last (which is a repeat of the first hour of the following water year) for each WY

    files = sorted(glob.glob('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_d02_' + str(year)+ '-' + str(m) + '-*'))[:-1]
    #files = sorted(glob.glob('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_*'))[:-1]
    print('/home/kmurenbeeld/scratch/wrf_out_to_subset/wy_' + str(year) + '/d02/wrfout_d02_' + str(year) + '-'+ str(m) + '-*')
    print("Total number of files: {}".format(len(files)))

    # open multi-file dataset (this function accepts unix wildcards)
    d = xr.open_mfdataset(files, drop_variables=var_list, concat_dim='Time', parallel=True).chunk({'Time':240})
    print(d)
    # Swap time and XTIME
    d = d.swap_dims({'Time':'XTIME'})

    # Hourly precip depth to rate
    d['PRCP'] = d['RAINNC']
    d['PRCP'].data = calc_precip(d['RAINNC'],d['I_RAINNC'])
    d['PRCP'] = d['PRCP']/3600

    # Calculate wind speed from the U and V components
    d['WIND'] = np.sqrt(d['U10']**2 + d['V10']**2)

    # Get 3-hourly data for the desired variables
    new_array = d[['T2','Q2','PSFC','SWNORM','GLW','PRCP','WIND']].resample(XTIME = '1H').mean(dim = 'XTIME') # create daily means of few variables
	new_array = new_array.rename({'T2':'TBOT','Q2':'QBOT','PSFC':'PSRF','SWNORM':'FSDS','GLW':'FLDS'}) # rename to be consistent with CLM(FATES)

	    # Adjust the meta data
	    new_array['TBOT'].attrs = [('description','1-HOURLY MEAN GRID SCALE TEMPERATUTE'), ('units','K')]
	    new_array['QBOT'].attrs = [('description','1-HOURLY MEAN GRID SCALE SPECIFIC HUMIDITY'), ('units','kg/kg')]
	    new_array['PSRF'].attrs = [('description','1-HOURLY MEAN SURFACE PRESSURE'), ('units','Pa')]
	    new_array['FSDS'].attrs = [('description','1-HOURLY MEAN NORMAL SHORT WAVE FLUX AT GROUND SURFACE (SLOPE-DEPENDENT)'), ('units','W m^2')]
	    new_array['FLDS'].attrs = [('description','1-HOURLY MEAN LONG WAVE FLUX AT GROUNF SURFACE'),('units','Wm^2')]
	    new_array['WIND'].attrs = [('description','1-HOURLY MAXIMUM GRID SCALE TEMPERATURE'), ('units','K')]
	    new_array['PRCP'].attrs = [('description','1-HOURLY MEAN PRECIPITATION'),('units','mm/s')]

	    # Write new netcdf file
	    #new_array.to_netcdf("/mnt/selway/data/data_02/charlie/subsets/test/forLejo/Biome-BGC-WY-" + str(year) + ".nc")
	    #new_array['PRCP'].to_netcdf(str(USER)+"@data-access.ucar.edu:/glade/scratch/katiem/clmforc/clmforc.WRF30d02.c2018.1kmx1km.Precip." + str(year) + "-" + str(month) + ".nc")
	    print("Writing output to disk")
	    new_array['PRCP'].to_netcdf("/home/kmurenbeeld/scratch/wrf_to_clm/clmforc_WRF30d02_c2018_1kmx1km_" + str(year) + "_" + str(m) + "_Precip.nc")
	    new_array['FSDS'].to_netcdf("/home/kmurenbeeld/scratch/wrf_to_clm/clmforc_WRF30d02_c2018_1kmx1km_" + str(year) + "_" + str(m) +  "_Solar.nc")
	  new_array[['TBOT','QBOT','PSRF','WIND','FLDS']].to_netcdf("/home/kmurenbeeld/scratch/wrf_to_clm/clmforc_WRF30d02_c2018_1kmx1km_" + str(year) + "_" + str(m) +  "_TPQWL.nc")

      t1 = time.time()
print("Total time to create this subset was:", t1 - t0, "seconds.")
