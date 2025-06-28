import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib.ticker as mticker
from matplotlib import gridspec
import pandas as pd

plt.close('all')

adt = xr.open_dataset('Ariane_workplace/currents_data/dt_global_allsat_phy_l4_20230428_20240501.nc')

lon_min_adt, lon_max_adt = 4,6
lat_min_adt, lat_max_adt = 40, 42
area_adt = [lon_min_adt, lon_max_adt, lat_min_adt, lat_max_adt]
adt = adt.sel(longitude=slice(lon_min_adt, lon_max_adt), latitude=slice(lat_min_adt, lat_max_adt))

adt_lon = adt['longitude'].values
adt_lon = adt['latitude'].values
adt_values = adt['adt'].squeeze().values

tsg_transect  = pd.read_csv('inputs/tsg_transect_AFB.csv')
tsg_station  = pd.read_csv('inputs/tsg_stations_AFB.csv')
A2 = tsg_station[tsg_station['Region'] == 'A2']
B2 = tsg_station[tsg_station['Region'] == 'B2']
F2 = tsg_station[tsg_station['Region'] == 'F2']
region_colors = {'A': 'lightblue', 'F': 'lightcoral', 'B': 'lightgreen'}
tsg_transect['color'] = tsg_transect['Region'].map(region_colors)
region_colors = {'A2': 'lightblue', 'F2': 'lightcoral', 'B2': 'lightgreen'}
tsg_station['color'] = tsg_station['Region'].map(region_colors)


def configure_map(ax, area):
    ax.set_extent(area, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, edgecolor='lightsteelblue')
    ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='lightsteelblue')
    ax.add_feature(cfeature.LAND, edgecolor='lightsteelblue', facecolor='ghostwhite')
    gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0, linewidth=0.5)
    gl.top_labels = False
    gl.bottom_labels = True
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(range(int(area[0]), int(area[1]) + 1))
    gl.ylocator = mticker.FixedLocator(range(int(area[2]), int(area[3]) + 1))
    gl.xformatter = mticker.FuncFormatter(lambda x, _: f"{abs(int(x))}°{'E' if x >= 0 else 'W'}")
    gl.yformatter = mticker.FuncFormatter(lambda y, _: f"{abs(int(y))}°{'N' if y >= 0 else 'S'}")


fig = plt.figure(figsize=(14, 6))
ax1 = fig.add_subplot(projection=ccrs.PlateCarree())
configure_map(ax1, area_adt)
img1 = plt.imshow(
    adt_values,
    extent=[lon_min_adt, lon_max_adt, lat_min_adt, lat_max_adt],
    origin='lower',
    cmap=cmocean.cm.topo,
    transform=ccrs.PlateCarree(),
    interpolation='nearest'
)

cbar = fig.colorbar(img1, ax=ax1, orientation='vertical', shrink=0.7, pad=0.05)
cbar.set_label('ADT (m)')

scatter = plt.scatter(tsg_transect['Longitude'], tsg_transect['Latitude'], c=tsg_transect['color'], s=8)
plt.scatter(A2['Longitude'], A2['Latitude'], color='lightblue', edgecolor='lightblue', s=8)
plt.scatter(F2['Longitude'], F2['Latitude'], color='lightcoral', edgecolor='lightcoral', s=8)
plt.scatter(B2['Longitude'], B2['Latitude'], color='lightgreen', edgecolor='lightgreen', s=8)

plt.scatter(np.mean(A2['Longitude']), np.mean(A2['Latitude']), color='lightblue', edgecolor='blue', marker='*', s=200, label="St.A2")
plt.scatter(np.mean(F2['Longitude']), np.mean(F2['Latitude']), color='lightcoral', edgecolor='red', marker='*', s=200, label="St.F2")
plt.scatter(np.mean(B2['Longitude']), np.mean(B2['Latitude']), color='lightgreen', edgecolor='green', marker='*', s=200, label="St.B2")

plt.show()

plt.savefig("outputs/adt_AB.tif", format="tiff", dpi=600)