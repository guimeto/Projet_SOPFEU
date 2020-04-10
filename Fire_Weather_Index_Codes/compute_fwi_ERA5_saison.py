# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:51:17 2019

@author: guillaume
"""
from fwi_compute_index import FWICLASS
from netCDF4 import Dataset
import numpy as np
import xarray as xr 
import pandas as pd 
import netCDF4
import gc 


gc.collect()


year = 2018 
multi_pr='J:/REANALYSES/ERA5/Prec_daily/Daily_Total_PR_17h_QC_'+str(year)+'*.nc'
multi_air='J:/REANALYSES/ERA5/HR/ERA5_QC_TC_HR_WIND_'+str(year)+'*.nc'

prec = xr.open_mfdataset(multi_pr).tp
hum = xr.open_mfdataset(multi_air).Humidity
wind = xr.open_mfdataset(multi_air).wind 
temp = xr.open_mfdataset(multi_air).Tc 
DS = xr.merge([prec,hum, wind,temp])
DS.to_netcdf('./tmp.nc')
del [prec, wind, hum, temp, DS] 

ds1 = Dataset('./tmp.nc')

prec = ds1.variables['tp'][:]
hum = ds1.variables['Humidity'][:]
wind = ds1.variables['wind'][:]
temp = ds1.variables['Tc'][:]
lats=ds1.variables['latitude'][:]
lons=ds1.variables['longitude'][:]

IND = np.zeros((len(prec),len(prec[0]),len(prec[0][0])),dtype=float)

dtime = netCDF4.num2date(ds1.variables['time'][:],ds1.variables['time'].units)

for ni in range(0, len(prec[0])):
    for nj in range(0, len(prec[0][0])):
        ffmc0 = 85.0
        dmc0 = 6.0
        dc0 = 15.0        
        for nt in range(0, len(prec)):
            fwisystem = FWICLASS(temp[nt,ni,nj],hum[nt,ni,nj],wind[nt,ni,nj],prec[nt,ni,nj])
            ffmc = fwisystem.FFMCcalc(ffmc0)
            dmc = fwisystem.DMCcalc(dmc0,dtime[nt].month)
            dc = fwisystem.DCcalc(dc0,dtime[nt].month)
            isi = fwisystem.ISIcalc(ffmc)
            bui = fwisystem.BUIcalc(dmc,dc)
            fwi = fwisystem.FWIcalc(isi,bui)
            ffmc0 = ffmc
            dmc0 = dmc
            dc0 = dc
            IND[nt,ni,nj]= fwi
 

data_set = xr.Dataset( coords={'lon': ([ 'lon'], ds1.variables['longitude'][:]),
                                 'lat': (['lat',], ds1.variables['latitude'][:]),
                                 'time': pd.date_range(str(year)+'-'+str(dtime[0].month)+'-01', periods=len(prec), freq='D')})
    
data_set["FWI"] = (['time','lat', 'lon'],  IND)
data_set = data_set.assign_coords(lon=(((data_set.lon + 180) % 360) - 180)).sortby('lon')

## load mask quebec
mask = xr.open_mfdataset('K:/PROJETS/PROJET_FIRE_INDEX/scripts/mask_QUEBEC.nc')
data_set=data_set.where(mask.region == 1)
del [ prec, hum, wind, temp, lats, lons] 

data_set.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_FWI_QC_'+str(year)+'_from_'+str(dtime[0].month)+'_to_'+str(dtime[nt].month)+'.nc')

# compute of severity rating SR = 0.0272 * FWI^1.77
data_set = 0.0272 * data_set 
SR = 0.0272 * np.power(data_set, 1.77)
SR = SR.resample(time = '1M').sum() 

SR.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_SR_QC_'+str(year)+'_from_'+str(dtime[0].month)+'_to_'+str(dtime[nt].month)+'.nc')

del [data_set, SR] 
 
                
