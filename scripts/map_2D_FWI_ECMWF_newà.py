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

dset=xr.open_mfdataset('K:/PROJETS/PROJET_FIRE_INDEX/FWI_ERA5_ECMWF/ECMWF_FWI_FWI_20180801_1200_hr_v3.1_con.nc')
#dtime = netCDF4.num2date(dset.variables['time'][:],dset.variables['time'].units)
#dset = dset.assign_coords(longitude=(((dset.longitude + 180) % 360) - 180)).sortby('longitude')

lat_bnd = [62, 43]
lon_bnd = [276, 306]
dset = dset.sel(longitude=slice(*lon_bnd), latitude=slice(*lat_bnd))
    

# select a day to display
year = 2018 
month_to_find = 8
day_to_find = 1 


   
## Lecture du fichier 
var=dset.fwi[:]
lon=dset.longitude[:]
lat=dset.latitude[:]
time = dset.time



fig = plt.figure(figsize=(28,16))
ax = plt.subplot(111, projection=ccrs.LambertConformal())
ax.set_extent([-90,-50,40,55])
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

## Choisissons une colormap
cmap0=plt.cm.jet
cmap0.set_under('darkblue') ## on met en blanc les valeurs inferieures au min de clev
cmap0.set_over('darkred') ## bleu fonce pour les valeurs extremes de pluie

mm = ax.pcolormesh(lon,\
                   lat,\
                   var.values,\
                   vmin=0,\
                   vmax=20, \
                   transform=ccrs.PlateCarree(),\
                   cmap=cmap0 )


  #  ax.stock_img();
mtl_lon, mtl_lat = -73.5, 45.5

# Define gridline locations and draw the lines using cartopy's built-in gridliner:
xticks = np.arange(-150.0,-40.0,20)
yticks =np.arange(10,80,10)

fig.canvas.draw()

# Standard 6,000 km scale bar.
scale_bar(ax, (0.70, 0.1), 100 ,plot_kwargs = dict(linestyle='dashed', color='black'))
cbar = plt.colorbar(mm,  shrink=0.75, drawedges='True', ticks=np.arange(0, 20.1, 1.), extend='both')
cbar.ax.tick_params(labelsize=20) 

string_title=u'Fire Weather Index: '
plt.title(string_title, size='xx-large')
plt.savefig('./FWI_ECMWF_new.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()






