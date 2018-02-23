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

fname_wrf = os.path.join('/', 'Users', 'tim', 'work', 'Data', 'WRF_Driver',
                         'met_em.d01.2009-06-01_00:00:00.nc')
lonwrf_ctr, latwrf_ctr = wrf_grid.get_wrf_latlon(fname_wrf, 'CLONG', 'CLAT')

# def in_wrf_domain(lonwrf, latwrf):
    # """test whether points in WRF domain are within the Redwoods range
    # """
fname = os.path.join('/', 'Users', 'tim', 'work', 'Data',
                     'Redwoods',
                     'Redwood_SequoiaSempervirens_extentNorthAmerica',
                     'data', 'commondata', 'data0', 'sequsemp.shp')

rwProj = Proj(init='epsg:3857')
latlonProj = Proj(init='epsg:4326')
xwrf_ctr, ywrf_ctr = transform(latlonProj, rwProj,
                               lonwrf_ctr.data, latwrf_ctr.data)

# bbrwsp = Point(-122.2464, 37.1737)  # Big Basin Redwoods State Park
# x2, y2 = transform(latlonProj, rwProj, -122.2464, 37.1737)
# bbrwsp = Point(x2, y2)  # Big Basin Redwoods State Park

rw_polys = fiona.open(fname)  # Redwood range polygons
rw_shapes = [shape(this_poly['geometry']) for this_poly in rw_polys]

has_redwoods = np.zeros(latwrf_ctr.shape, dtype=bool)
for i in range(latwrf_ctr.shape[0]):
    for j in range(latwrf_ctr.shape[1]):
        p = Point(xwrf_ctr[i, j], ywrf_ctr[i, j])
        has_redwoods[i, j] = np.any([p.within(this_shape) for
                                     this_shape in rw_shapes])

(lon_LL, lat_LL,
 lon_UL, lat_UL,
 lon_UR, lat_UR,
 lon_LR, lat_LR) = wrf_grid.get_cell_corner_coords(fname_wrf)
x_LL, y_LL = transform(latlonProj, rwProj, lon_LL, lat_LL)

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
has_redwoods = ma.masked_where(has_redwoods == False, has_redwoods)

ax.pcolormesh(x_LL, y_LL,
              has_redwoods,
              transform=ccrs.Mercator.GOOGLE)
ax.scatter(xwrf_ctr, ywrf_ctr,
           color='black', marker='x',
           transform=ccrs.Mercator.GOOGLE,
           label='center')
rw_shapes_c = list(shpreader.Reader(fname).geometries())
ax.add_geometries(geoms=rw_shapes_c, crs=ccrs.Mercator.GOOGLE,
                  edgecolor='blue', facecolor='gray', alpha=0.5)

ax.legend()
plt.show()
