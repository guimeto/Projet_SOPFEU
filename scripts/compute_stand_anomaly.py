# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:23:46 2020

@author: guillaume
"""

import xarray as xr 
####https://uoftcoders.github.io/studyGroup/lessons/python/cartography/lesson/

# lecture de la serie complete
file = 'K:/PROJETS/PROJET_FIRE_INDEX/ERA5_Onset_Netcdf/ERA5_SNOW_Onset_'
multi_file = [f'{file}{year}_from_3_to_8.nc' for year in range(1979,2020,1)]
ds_all = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')

data_all = ds_all.variables['Onset'][:].squeeze()


# lecture de la serie pour calculer la clim et la deviation
file = 'K:/PROJETS/PROJET_FIRE_INDEX/ERA5_Onset_Netcdf/ERA5_SNOW_Onset_'
multi_file = [f'{file}{year}_from_3_to_8.nc' for year in range(1990,2020,1)]

ds_clim = xr.concat([xr.open_dataset(f) for f in multi_file], 'time')

data_clim = ds_clim.variables['Onset'][:].mean("time")
data_std = ds_clim.variables['Onset'][:].std("time")

stand_anomalies = xr.apply_ufunc(
    lambda x, m, s: (x - m) / s,
    ds_all.groupby("time"),
    data_clim,
    data_std,
)

stand_anomalies.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_Onset_Netcdf/ERA5_Onset_Anomaly_Stand_1979-2019_vs_1990-2019.nc')








