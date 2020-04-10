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

dset=Dataset('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_FWI_QC2_2018_from_1_to_12.nc')
dtime = netCDF4.num2date(dset.variables['time'][:],dset.variables['time'].units)

# select a day to display
year = 2018 
month_to_find = 8
day_to_find = 1 

for index, item in enumerate(dtime):
    if (dtime[index].month == month_to_find and  dtime[index].day == day_to_find):
        index_to_find = index   
print('On va tracer la journée du: ' + dtime[index_to_find].strftime("%Y-%m-%d") )

## Lecture du fichier 
var=dset.variables['FWI'][index_to_find][:]
lon=dset.variables['lon'][:]
lat=dset.variables['lat'][:]
time = dset.variables['time']

data_set = xr.Dataset( coords={'lon': ([ 'x'], lon),
                                 'lat': (['y',], lat),
                                 'time': dtime[index_to_find].strftime("%Y-%m-%d")})
    

    
data_set.to_netcdf('K:/PROJETS/PROJET_FIRE_INDEX/ERA5_FWI_Netcdf/ERA5_FWI_QC_'+dtime[index_to_find].strftime("%Y-%m-%d")+'.nc')

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
                   var,\
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

string_title=u'Fire Weather Index: '+dtime[index_to_find].strftime("%Y-%m-%d") 
plt.title(string_title, size='xx-large')
plt.savefig('./FWI_bis2_'+dtime[index_to_find].strftime("%Y-%m-%d")+'.png', bbox_inches='tight', pad_inches=0.1)
plt.show()  
plt.close()






