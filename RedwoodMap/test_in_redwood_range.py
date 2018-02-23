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

def get_wrf_latlon(wrf_file, lonvar='XLONG', latvar='XLAT'):
    """read latitude and longitude from WRF input or output
    """
    nc = netCDF4.Dataset(wrf_file)
    lat = getvar(nc, latvar)
    lon = getvar(nc, lonvar)
    nc.close()
    return(lon, lat)

def half_cell_shift(x, y):
    dx = np.diff(x, axis=1)
    dy = np.diff(y, axis=0)
    xshift = x.copy()
    xshift[:, :-1] -= (dx / 2.0)
    yshift = y.copy()
    yshift[:-1, :] -= (dy / 2.0)
    return(xshift, yshift)

fname_wrf = os.path.join('/', 'Users', 'tim', 'work', 'Data', 'WRF_Driver',
                         'met_em.d01.2009-06-01_00:00:00.nc')
lonwrf_crnr, latwrf_crnr = get_wrf_latlon(fname_wrf, 'XLONG_C', 'XLAT_C')
lonwrf_ctr, latwrf_ctr = get_wrf_latlon(fname_wrf, 'CLONG', 'CLAT')

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
xwrf_crnr, ywrf_crnr = transform(latlonProj, rwProj,
                                 lonwrf_crnr.data, latwrf_crnr.data)

rw_polys = fiona.open(fname)  # Redwood range polygons
rw_shapes = [shape(this_poly['geometry']) for this_poly in rw_polys]

# wrf_cell_centers = Point(xwrf, ywrf)

# for i, this_shape in enumerate(rw_shapes):
#     in_poly = bbrwsp.within(this_shape)
#     print("{:02d} {}".format(i, in_poly))





# bbrwsp = Point(-122.2464, 37.1737)  # Big Basin Redwoods State Park
# x2, y2 = transform(latlonProj, rwProj, -122.2464, 37.1737)
# bbrwsp = Point(x2, y2)  # Big Basin Redwoods State Park

has_redwoods = np.zeros(latwrf_ctr.shape, dtype=bool)
for i in range(latwrf_ctr.shape[0]):
    for j in range(latwrf_ctr.shape[1]):
        p = Point(xwrf_ctr[i, j], ywrf_ctr[i, j])
        has_redwoods[i, j] = np.any([p.within(this_shape) for
                                     this_shape in rw_shapes])


# def myPoint(x, y):
#     return(Point(x, y))
# pointvec = np.vectorize(myPoint)
# p2 = np.empty(latwrf.shape, dtype=object)
# p2 = pointvec(xwrf, ywrf)

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
ax.set_extent((lonwrf_ctr.min(), lonwrf_ctr.max(),
               latwrf_ctr.min(), latwrf_ctr.max()))
has_redwoods = ma.masked_where(has_redwoods == False, has_redwoods)
ax.pcolormesh(xwrf_crnr, ywrf_crnr, has_redwoods,
              transform=ccrs.Mercator.GOOGLE)
ax.scatter(xwrf_ctr, ywrf_ctr,
           color='black', marker='x',
           transform=ccrs.Mercator.GOOGLE)
rw_shapes_c = list(shpreader.Reader(fname).geometries())
ax.add_geometries(geoms=rw_shapes_c, crs=ccrs.Mercator.GOOGLE,
                  edgecolor='blue', facecolor='gray', alpha=0.5)

plt.show()
