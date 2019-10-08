# coding: utf-8

import numpy as np
import xarray as xr
import pandas as pd
import regrid
import dask.array as da
import dask

"""Regrid 0.5x0.625 MERRA-2 soil wetness to 2x2.5 resolution"""
# Takes ~2 mins for a month

# ## Set up coordinates for regridding
D2R = 3.141592658979323 / 180.0

DI = 2.5
DJ = 2.
nI = 144
nJ = 91
lon_b = np.array([-180 - DI/2 + ( DI * I ) for I in np.arange(nI+1)])
lat_b = np.array([-90 - DJ/2 + ( DJ * J ) for J in np.arange(nJ+1)])
lat_b[[0,-1]] = [-90,90] # reset poles
lon = np.array([-180 + ( DI * I ) for I in np.arange(nI)])
lat = np.array([-90 + ( DJ * J ) for J in np.arange(nJ)])
lat[[0,-1]] = [-89.75,89.75]
sine_lat = np.sin( lat_b * D2R )

ds_low = xr.Dataset(
    None, 
    coords = {
        'lat': lat,
        'lon': lon,
        'lat_b': lat_b,
        'lon_b': lon_b,
        'sin_lat_b': ('lat_b', sine_lat),
    }
)

DI = 0.625
DJ = 0.5
nI = 576
nJ = 361
lon_b = np.array([-180 - DI/2 + ( DI * I ) for I in np.arange(nI+1)])
lat_b = np.array([-90 - DJ/2 + ( DJ * J ) for J in np.arange(nJ+1)])
lat_b[[0,-1]] = [-90,90] # reset poles
lon = np.array([-180 + ( DI * I ) for I in np.arange(nI)])
lat = np.array([-90 + ( DJ * J ) for J in np.arange(nJ)])
lat[[0,-1]] = [-89.875,89.875]
sine_lat = np.sin( lat_b * D2R )

ds_high = xr.Dataset(
    None, 
    coords = {
        'lat': lat,
        'lon': lon,
        'lat_b': lat_b,
        'lon_b': lon_b,
        'sin_lat_b': ('lat_b', sine_lat),
    }
)


def h2l_wrapper(high_res):
    """Regrids from 0.5x0.625 to 2x2.5, mimicking how the GCST regrids MERRA-2 data"""
    values = regrid.regridmodule.map_a2a( 
        ds_high.lon_b.values, ds_high.sin_lat_b.values, 
        high_res.transpose(), 
        ds_low.lon_b.values, ds_low.sin_lat_b.values,           
        0, 0,
        576, 361, 144, 91, 
    )
    return values.transpose()

# ## Read data from file

year = 2012
time = pd.date_range(start='{:4d}0101 00:30:00'.format(year), end='{:4d}1231 23:30:00'.format(year), freq='H')
# for month in range(1,13):
month = 1
dates = time[time.month == month][::24].strftime('%Y%m%d')

# sample data
low_res_sample = xr.open_dataset('/Users/tgabi/Desktop/Data/MERRA2/MERRA2.20170101.A1.2x25.nc4')[['GWETROOT','GWETTOP']].isel(time=0)

# high resolution data to be regridded
files = ['/Users/tgabi/Desktop/Data/MERRA2_soil_moisture/MERRA2_400.tavg1_2d_lnd_Nx.{}.nc4.nc'.format(date) for date in dates]
# files = ['/Users/tgabi/Desktop/Data/MERRA2_soil_moisture/MERRA2_400.tavg1_2d_lnd_Nx.{}.nc4.nc'.format(dates[0])]
list_objs = [xr.open_dataset(fname) for fname in files]
vars_to_include = ['GWETROOT','GWETTOP']

# weighting applied in preprocessing 
weight = (xr.open_dataset('/Users/tgabi/Desktop/Data/MERRA2_101.const_2d_asm_Nx.00000000.nc4')['FRLAND'].isel(time=0) * 
          xr.open_dataset('/Users/tgabi/Desktop/Data/MERRA2_100.const_2d_lnd_Nx.00000000.nc4')['poros'].isel(time=0)
         )
weight = weight.where(np.isfinite(weight),0).values

print('regridding...')

# pre-processing: set value=0 over grids with zero land fraction
tmp_arrays = [da.concatenate([da.from_array(obj[var].where(weight!=0,0), chunks='auto') for obj in list_objs], axis=0) for var in vars_to_include]
# pre-processing: multiply by porosity and land fraction
denominator = h2l_wrapper(weight)
# regridding with weights applied
arrays = [da.stack([h2l_wrapper(array[i,:,:]*weight) / denominator for i in range(array.shape[0])], axis=0) for array in tmp_arrays]
computed = dask.compute(*arrays)

# save to dataset
h2l = xr.Dataset(dict([(var, (('time','lat','lon'), array)) for var,array in zip(vars_to_include, computed)]), 
                    coords = {'time':time[time.month == month],
                    # coords = {'time':time[time.month == month][:24],
                              'lat': low_res_sample.lat,
                              'lon': low_res_sample.lon,
                             }
                   )

# post-processing: set nans to 0. They arise from 0/0 when land fraction is 0
h2l['GWETROOT'] = h2l['GWETROOT'].where(denominator!=0, 0)
h2l['GWETTOP'] = h2l['GWETTOP'].where(denominator!=0, 1)
# set maximum value to 1
h2l['GWETROOT'] = h2l['GWETROOT'].where(h2l['GWETROOT'] < 1, 1)
h2l['GWETTOP'] = h2l['GWETTOP'].where(h2l['GWETTOP'] < 1, 1)

# add attributes
h2l['GWETROOT'].attrs = low_res_sample['GWETROOT'].attrs
h2l['GWETTOP'].attrs = low_res_sample['GWETTOP'].attrs

# save to netcdf
h2l.to_netcdf('/Users/tgabi/Desktop/Data/MERRA2_soil_moisture/MERRA2_modified_wetness.{:4d}{:02d}.2x25.nc4'.format(year, month))
