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

PATH_TO_SHAPEFILE = './gpr_000b11a_f/gpr_000b11a_f.shp'
prov = gpd.read_file(PATH_TO_SHAPEFILE)
prov.head()

data = 'ERA5_FWI_QC_201707.nc'
ds = xr.open_mfdataset(data, chunks = {'time': 10})
ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
ds.lon.values

ID_PROV = 1
print(prov.PRFNOM[ID_PROV])

prov_mask_poly = regionmask.Regions_cls(name = 'PRFNOM', numbers = list(range(0,3)), names = list(prov.PRFNOM), abbrevs = list(prov.PRFNOM), outlines = list(prov.geometry.values[i] for i in range(0,3)))
prov_mask_poly


mask = prov_mask_poly.mask(ds.isel(time = 0), lat_name='lat', lon_name='lon')

mask.to_netcdf('./mask_QUEBEC.nc') 
