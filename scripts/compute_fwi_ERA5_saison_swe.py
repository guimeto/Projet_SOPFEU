# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:51:17 2019

@author: guillaume

Calcul du FWI avec condition de neige en surface
L’initialisation doit se faire en début de saison (habituellement 3-2 jours
 après la fonte de la neige) et les indices se calcul de manière continue 
par la suite pour toute la saison. 

"""
from fwi_compute_index import FWICLASS
from netCDF4 import Dataset
import numpy as np
import xarray as xr 
import pandas as pd 
import netCDF4
import gc 


gc.collect()

yi = 2018
yf = 2018
#########################################################

for year in range(yi,yf+1):
    print('Calcul pour: '+str(year))
    # ouverture des fichiers pour une année 
    # precipitation cumulée sur 24h
    file = 'J:/REANALYSES/ERA5/Prec_daily/Daily_Total_PR_17h_QC_'+str(year)
    multi_pr = [f'{file}{month}.nc' for month in ['04','05','06','07','08']]
    # variables humidité, temperature et vent à 17h UTC
    file = 'J:/REANALYSES/ERA5/HR/ERA5_QC_TC_HR2_WIND_'+str(year)
    multi_air = [f'{file}{month}_17h.nc' for month in ['04','05','06','07','08']]
    # épaisseur de neige quotidienne 
    multi_snow = ['J:/REANALYSES/ERA5/SWE/ERA5_'+str(year)+month+'/ERA5_SWE_'+str(year)+month+'_daymean.nc4' for month in ['04','05','06','07','08']]
    swe = xr.open_mfdataset(multi_snow)
    
    # on coupe les données de neige au dessus du Quebec
    lat_bnd = [62, 43]
    lon_bnd = [276, 306]
    sd = swe.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd)).sd
         
    prec = xr.open_mfdataset(multi_pr).tp
    hum = xr.open_mfdataset(multi_air).Humidity
    wind = xr.open_mfdataset(multi_air).wind 
    temp = xr.open_mfdataset(multi_air).Tc 
    DS = xr.merge([prec,hum, wind,temp])
    
    DS.to_netcdf('./tmp.nc')
    del [prec, wind, hum, temp, DS] 
    ds1 = Dataset('./tmp.nc')
    
    # lecture des champs pour le script de calcul du FWI
    
    prec = ds1.variables['tp'][:]
    hum = ds1.variables['Humidity'][:]
    wind = ds1.variables['wind'][:]
    temp = ds1.variables['Tc'][:]
    snow = sd.values
    
    lats=ds1.variables['latitude'][:]
    lons=ds1.variables['longitude'][:]
    
    # declaration du champs de sortie:  nlats X nlons X durée saison
    IND = np.zeros((len(prec),len(prec[0]),len(prec[0][0])),dtype=float)
    dtime = netCDF4.num2date(ds1.variables['time'][:],ds1.variables['time'].units)
    
    # travail en chaque point de grille
    for ni in range(0, len(prec[0])):
        for nj in range(0, len(prec[0][0])):
            
            # detection du premier jour sans neige au sol en un point donné
            ini_t = list(snow[:,ni,nj]).index(0)
            
            # si plus de neige en surface au temps t, on initialise les calculs du FWI à t+3
            ffmc0 = 85.0
            dmc0 = 6.0
            dc0 = 15.0 
            
            for nt in range(ini_t+3, len(prec)):
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
    
    # sauvegarde du champs en xarray
    data_set = xr.Dataset( coords={'lon': ([ 'lon'], ds1.variables['longitude'][:]),
                                     'lat': (['lat',], ds1.variables['latitude'][:]),
                                     'time': pd.date_range(str(year)+'-'+str(dtime[0].month)+'-01', periods=len(prec), freq='D')})
        
    data_set["FWI"] = (['time','lat', 'lon'],  IND)
    data_set = data_set.assign_coords(lon=(((data_set.lon + 180) % 360) - 180)).sortby('lon')
    
    ## load mask quebec
    mask = xr.open_mfdataset('K:/PROJETS/PROJET_FIRE_INDEX/scripts/mask_QUEBEC.nc')
    data_set=data_set.where(mask.region == 1)
    del [ prec, hum, wind, temp, lats, lons] 
    
    # ecriture du Netcdf
    data_set.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_FWI_QC_SNOW1old_'+str(year)+'_from_'+str(dtime[0].month)+'_to_'+str(dtime[nt].month)+'.nc')
    
    # compute of severity rating SR = 0.0272 * FWI^1.77 
    SR = 0.0272 * np.power(data_set, 1.77)
    # pour cumuler cet indeice de sévérité sur une période mensuelle ou saisonnière:
    #SR = SR.resample(time = '1M').sum() 
    # pour cumuler sur l'ensemble de la saison 
    SR = SR.sum(dim='time')

    SR.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_SR_QC_SNOW1old_'+str(year)+'_from_'+str(dtime[0].month)+'_to_'+str(dtime[nt].month)+'.nc')
  
    ds1.close()
    del [data_set, SR] 

              
