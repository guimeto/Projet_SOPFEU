# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:51:17 2019

@author: guillaume

Calcul du FWI avec condition de neige en surface
L’initialisation doit se faire en début de saison (habituellement 3-2 jours
 après la fonte de la neige) et les indices se calcul de manière continue 
par la suite pour toute la saison. 

"""
from netCDF4 import Dataset
import numpy as np
import xarray as xr 
import pandas as pd 
import netCDF4
import gc 
import math


"""     Define Class FWICLASS first   """

class FWICLASS:
    def __init__(self,temp,rhum,wind,prcp):
        self.h=rhum
        self.t=temp
        self.w=wind
        self.p=prcp
        
    def FFMCcalc(self,ffmc0): 
        global m
        mo = (147.2*(101.0 - ffmc0))/(59.5 + ffmc0) #*Eq.1*#
        if (self.p > 0.5):
                rf = self.p - 0.5                   #*Eq.2*#
                if(mo > 150.0):
                        mo=(mo+42.5*rf*math.exp(-100.0/(251.0-mo))*(1.0 - math.exp(-6.93/rf))) \
                 + (.0015*(mo - 150.0)**2)*math.sqrt(rf) #*Eq.3b*#
                elif mo<=150.0:
                        mo=mo+42.5*rf*math.exp(-100.0/(251.0-mo))*(1.0 - math.exp(-6.93/rf))
                                                                         #*Eq.3a*#
                if(mo > 250.0):
                    mo=250.0 

        ed=.942*(self.h**.679) + (11.0*math.exp((self.h-100.0)/10.0))+0.18*(21.1-self.t) \
                *(1.0 - 1.0/math.exp(.1150 * self.h))    #*Eq.4*#

        if(mo < ed):
                ew = .618*(self.h**.753) + (10.0*math.exp((self.h-100.0)/10.0)) \
             + .18*(21.1-self.t)*(1.0 - 1.0/math.exp(.115 * self.h)) #*Eq.5*#
                if(mo <= ew):
                        kl = .424*(1.0-((100.0-self.h)/100.0)**1.7)+(.0694*math.sqrt(self.w)) \
                                *(1.0 - ((100.0 - self.h)/100.0)**8) #*Eq.7a*#
                        kw = kl * (.581 * math.exp(.0365 * self.t))  #*Eq.7b*#
                        m  = ew - (ew - mo)/10.0**kw                 #*Eq.9*#
                elif mo > ew:
                        m = mo
        elif(mo == ed):
                m = mo
        elif (mo > ed):
                kl =.424*(1.0-(self.h/100.0)**1.7)+(.0694*math.sqrt(self.w))*   \
                        (1.0-(self.h/100.0)**8)  #*Eq.6a*#
                kw = kl * (.581*math.exp(.0365*self.t))          #*Eq.6b*#
                m  = ed + (mo-ed)/10.0 ** kw                     #*Eq.8*#

        ffmc = (59.5 * (250.0 -m)) / (147.2 + m)                 #*Eq.10*#
        if (ffmc  > 101.0):
            ffmc = 101.0
        if (ffmc  <=  0.0):
            ffmc =   0.0
        return ffmc

    def DMCcalc(self,dmc0,mth):
        global b
        el=[6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,6.0]

        t = self.t
        if(t < -1.1):
            t = -1.1
        rk = 1.894*(t+1.1) * (100.0-self.h) * (el[mth-1]*0.0001) #*Eqs.16 and 17*#

        if self.p > 1.5:
            ra= self.p 
            rw = 0.92*ra - 1.27                              #*Eq.11*#
            wmi=20.0 + 280.0/math.exp(0.023*dmc0)            #*Eq.12*#
            if dmc0 <= 33.0:
                b = 100.0 /(0.5 + 0.3*dmc0)              #*Eq.13a*#
            elif dmc0 > 33.0:
                if dmc0 <= 65.0: 
                    b = 14.0 - 1.3*math.log(dmc0)    #*Eq.13b*#
            elif dmc0 > 65.0:
                    b =  6.2 * math.log(dmc0) - 17.2 #*Eq.13c*#
            wmr = wmi + (1000*rw) / (48.77+b*rw)             #*Eq.14*#
            pr = 43.43 * (5.6348 - math.log(wmr-20.0))       #*Eq.15*#
        elif self.p <= 1.5:
            pr = dmc0 
        if(pr<0.0):
            pr=0.0
        dmc = pr + rk
        if(dmc<=1.0):
            dmc=1.0
        return dmc

    def DCcalc(self,dc0,mth):
        global dc
        fl=[-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6]
        t = self.t

        if(t < -2.8):
            t = -2.8  
        pe = (  0.36*(t+2.8) + fl[mth-1] )/2    #*Eq.22*#
        if pe<=0.0:
            pe = 0.0

        if (self.p > 2.8):
            ra = self.p
            rw = 0.83*ra - 1.27             #*Eq.18*#
            smi=800.0 * math.exp(-dc0/400.0)#*Eq.19*#
            dr = dc0 - 400.0*math.log(  1.0+((3.937*rw)/smi) ) #*Eqs. 20 and 21*#
            
            if(dr > 0.0):
                dc = dr + pe
        elif self.p <= 2.8:
            dc = dc0 + pe
        return dc

    def ISIcalc(self,ffmc):
        mo = 147.2*(101.0-ffmc) / (59.5+ffmc)                         #*Eq.1*#
        ff = 19.115*math.exp(mo*-0.1386) * (1.0+(mo**5.31)/49300000.0)#*Eq.25*#
        isi = ff * math.exp(0.05039*self.w)                           #*Eq.26*#
        return isi


    def BUIcalc(self,dmc,dc):
        if dmc <= 0.4*dc:
            bui = (0.8*dc*dmc) / (dmc+0.4*dc)                     #*Eq.27a*#
        else:
            bui=dmc-(1.0-0.8*dc/(dmc+0.4*dc))*(0.92+(0.0114*dmc)**1.7)#*Eq.27b*#
        if bui <0.0:
            bui=0.0
        return bui


    def FWIcalc(self,isi,bui):
        if bui <= 80.0:
            bb = 0.1 * isi * (0.626*bui**0.809 + 2.0)                #*Eq.28a*#
        else:
            bb = 0.1*isi*(1000.0/(25. + 108.64/math.exp(0.023*bui))) #*Eq.28b*#
        if(bb <= 1.0):
             fwi = bb                                                 #*Eq.30b*#
        else:
            fwi = math.exp(2.72 *  (0.434*math.log(bb))**0.647)      #*Eq.30a*#

        return fwi


"""    End of class FWICLASS    """

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
    # humidity
    file = 'J:/REANALYSES/ERA5/HR/ERA5_QC_Humdidity_'+str(year)
    multi_hum = [f'{file}{month}_17h.nc' for month in ['04','05','06','07','08']]
    # épaisseur de neige quotidienne 
    multi_snow = ['J:/REANALYSES/ERA5/SWE/ERA5_'+str(year)+month+'/ERA5_SWE_'+str(year)+month+'_daymean.nc4' for month in ['04','05','06','07','08']]
    swe = xr.open_mfdataset(multi_snow)
    
    # on coupe les données de neige au dessus du Quebec
    lat_bnd = [62, 43]
    lon_bnd = [276, 306]
    sd = swe.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd)).sd
         
    prec = xr.open_mfdataset(multi_pr).tp
    hum = xr.open_mfdataset(multi_hum).Humidity
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
    data_set.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_FWI_QC_SNOW_New_'+str(year)+'_from_'+str(dtime[0].month)+'_to_'+str(dtime[nt].month)+'.nc')
    
    # compute of severity rating SR = 0.0272 * FWI^1.77 
    SR = 0.0272 * np.power(data_set, 1.77)
    # pour cumuler cet indeice de sévérité sur une période mensuelle ou saisonnière:
    #SR = SR.resample(time = '1M').sum() 
    # pour cumuler sur l'ensemble de la saison 
    SR = SR.sum(dim='time')
    SR.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_SR_QC_SNOW_New_'+str(year)+'_from_'+str(dtime[0].month)+'_to_'+str(dtime[nt].month)+'.nc')
    
    ds1.close()
    del [data_set, SR] 

              
