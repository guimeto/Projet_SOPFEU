# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:23:46 2020

@author: guillaume
"""
import netCDF4
from netCDF4 import Dataset
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from carto import scale_bar
from datetime import datetime
import xarray as xr 
####https://uoftcoders.github.io/studyGroup/lessons/python/cartography/lesson/

## Date à utiliser 

dset=Dataset('K:\PROJETS\PROJET_FIRE_INDEX\ERA5_Onset_Netcdf/ERA5_SNOW_Onset_2018_from_4_to_8.nc')

year = 2018
## Lecture du fichier 
var=dset.variables['Onset'][:]
lon=dset.variables['lon'][:]
lat=dset.variables['lat'][:]


fig = plt.figure(figsize=(28,16))
ax = plt.subplot(111, projection=ccrs.LambertConformal())
ax.set_extent([-90,-50,40,60])
   # ax.coastlines(resolution='110m');
ax.add_feature(cfeature.OCEAN.with_scale('50m'))      # couche ocean
ax.add_feature(cfeature.LAND.with_scale('50m'))       # couche land
ax.add_feature(cfeature.LAKES.with_scale('50m'))      # couche lac    
ax.add_feature(cfeature.BORDERS.with_scale('50m'))    # couche frontieres
ax.add_feature(cfeature.RIVERS.with_scale('50m'))     # couche rivières 
coast = cfeature.NaturalEarthFeature(category='physical', scale='10m',     # ajout de la couche cotière 
                        facecolor='none', name='coastline')
ax.add_feature(coast, edgecolor='black')

  
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')


mm = ax.contourf(lon,\
                   lat,\
                   var,\
                   vmin=0,\
                   vmax=100, \
                   transform=ccrs.PlateCarree(),\
                   levels=np.arange(0, 100, 10.),\
                   cmap=plt.cm.jet )

  #  ax.stock_img();
mtl_lon, mtl_lat = -73.5, 45.5

# Define gridline locations and draw the lines using cartopy's built-in gridliner:
xticks = np.arange(-150.0,-40.0,20)
yticks =np.arange(10,80,10)

fig.canvas.draw()

# Standard 6,000 km scale bar.
scale_bar(ax, (0.70, 0.1), 100 ,plot_kwargs = dict(linestyle='dashed', color='black'))
cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(0, 400.1, 20.), extend='both')
cbar.ax.tick_params(labelsize=20) 

string_title=u'FWI Onset: season '+ str(year) 
plt.title(string_title, size='xx-large')
plt.savefig('./figures/Onset_FWI_'+str(year)+'.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()







