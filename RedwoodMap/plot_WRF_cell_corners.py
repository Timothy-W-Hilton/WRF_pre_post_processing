"""script to plot WRF cell centers and corners on a map

If the map looks reasonable it demonstrates that
wrf_grid.get_cell_corner_coords works.
"""

import numpy as np
import numpy.ma as ma
import os.path
import fiona
from shapely.geometry import shape,mapping, Point, Polygon, MultiPolygon
from pyproj import Proj, transform
import netCDF4
from wrf import getvar

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt

import wrf_grid

# ==================================================

fname_wrf = os.path.join('/', 'Users', 'tim', 'work', 'Data', 'WRF_Driver',
                         'met_em.d01.2009-06-01_00:00:00.nc')

(lon_LL, lat_LL,
 lon_UL, lat_UL,
 lon_UR, lat_UR,
 lon_LR, lat_LR) = wrf_grid.get_cell_corner_coords(fname_wrf)

lon_c, lat_c = wrf_grid.get_wrf_latlon(fname_wrf,
                                       'CLONG', 'CLAT')

fig = plt.figure(figsize=(6, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution='10m')
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')
ax.add_feature(states_provinces, edgecolor='grey')
# ax.add_geometries(geoms=rw_shapes, crs=proj,
#                   edgecolor='blue', facecolor='blue')
ax.set_extent((lon_LL.min(), lon_UR.max(),
               lat_LL.min(), lat_UR.max()))

ax.scatter(lon_c, lat_c,
           color='black', marker='x',
           transform=ccrs.PlateCarree(),
           label='center')
# colors from colorbrewer2.org 'dark2' palette
ax.scatter(lon_LL, lat_LL,
           color='#1b9e77', marker='x',
           transform=ccrs.PlateCarree(),
           label='LL')
ax.scatter(lon_UL, lat_UL,
           color='#d95f02', marker='x',
           transform=ccrs.PlateCarree(),
           label='UL')
ax.scatter(lon_UR, lat_UR,
           color='#7570b3', marker='x',
           transform=ccrs.PlateCarree(),
           label='UR')
ax.scatter(lon_LR, lat_LR,
           color='#e7298a', marker='x',
           transform=ccrs.PlateCarree(),
           label='LR')

ax.legend()
plt.show()
