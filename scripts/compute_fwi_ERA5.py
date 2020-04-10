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



year = 2018 
mth = 8

m = f"{mth:02d}"

ds1 = Dataset("J:/REANALYSES/ERA5/Prec_daily/Daily_Total_PR_17h_QC_"+str(year)+str(m)+".nc")
ds2 = Dataset("J:/REANALYSES/ERA5/HR/ERA5_QC_TC_HR_WIND_"+str(year)+str(m)+"_17h.nc")


prec = ds1.variables['tp'][:]
hum = ds2.variables['Humidity'][:]
wind = ds2.variables['wind'][:]
temp = ds2.variables['Tc'][:]
lats=ds1.variables['latitude'][:]
lons=ds1.variables['longitude'][:]

IND = np.zeros((len(prec),len(prec[0]),len(prec[0][0])),dtype=float)


for ni in range(0, len(prec[0])):
    for nj in range(0, len(prec[0][0])):
        ffmc0 = 85.0
        dmc0 = 6.0
        dc0 = 15.0
        
        for nt in range(0, len(prec)):
            fwisystem = FWICLASS(temp[nt,ni,nj],hum[nt,ni,nj],wind[nt,ni,nj],prec[nt,ni,nj])
            ffmc = fwisystem.FFMCcalc(ffmc0)
            dmc = fwisystem.DMCcalc(dmc0,mth)
            dc = fwisystem.DCcalc(dc0,mth)
            isi = fwisystem.ISIcalc(ffmc)
            bui = fwisystem.BUIcalc(dmc,dc)
            fwi = fwisystem.FWIcalc(isi,bui)
            ffmc0 = ffmc
            dmc0 = dmc
            dc0 = dc
            IND[nt,ni,nj]= fwi
 

data_set = xr.Dataset( coords={'lon': ([ 'x'], ds1.variables['longitude'][:]),
                                 'lat': (['y',], ds1.variables['latitude'][:]),
                                 'time': pd.date_range(str(year)+'-'+str(m)+'-01', periods=31, freq='D')})
    
data_set["FWI"] = (['time','y', 'x'],  IND)
    
data_set.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_FWI_QC_'+str(year)+str(m)+'.nc')
 
                
