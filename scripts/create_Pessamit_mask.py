# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import xarray as xr 
import numpy as np
import regionmask
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt

PATH_TO_SHAPEFILE = './shapes_pessamit/Nitassinan_essipit.shp'
shapes = gpd.read_file(PATH_TO_SHAPEFILE)
shapes
shapes.loc[0, 'geometry']
shapes.plot()

#tmp = gpd.GeoDataFrame.from_file('nitassinan_l.shp')
tmpWGS84 = shapes.to_crs({'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'})
tmpWGS84.loc[0, 'geometry']
tmpWGS84.to_file('Nitassinan_essipit_WGS84.shp')
tmpWGS84.plot()


data = 'K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_FWI_QC_SNOW_New_2019_from_4_to_8.nc'
ds = xr.open_mfdataset(data, chunks = {'time': 10})
ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')


pessamit_mask_poly = regionmask.Regions_cls(name = 'DESCRIPTIO', numbers = list(range(0,1)), names = list(tmpWGS84.DESCRIPTIO), abbrevs = list(tmpWGS84.DESCRIPTIO), outlines = list(tmpWGS84.geometry.values[i] for i in range(0,1)))



mask = pessamit_mask_poly.mask(ds.isel(time = 0), lat_name='lat', lon_name='lon')

mask.to_netcdf('./mask_Pessamit.nc') 
