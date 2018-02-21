"""modified from https://gis.stackexchange.com/questions/259831/shapefiles-wont-display-with-cartopy
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt

fname = r'/Users/tim/work/Data/Redwoods/Redwood_SequoiaSempervirens_extentNorthAmerica/data/commondata/data0/sequsemp.shp'
rw_shapes = list(shpreader.Reader(fname).geometries())

# a number of places on the web
# (e.g. https://github.com/SciTools/cartopy/issues/924) suggest
# crs.epsg(3857) and crs.Mercator.GOOGLE are the same.  But the
# redwood range plots differently.  crs.Mercator.GOOGLE.proj4_params
# and crs.epsg(3857).proj4_params differ; that's another way to
# illustrate this.
#proj = ccrs.epsg(3857)
proj = ccrs.Mercator.GOOGLE

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution='10m')
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')
ax.add_feature(states_provinces, edgecolor='grey')
ax.add_geometries(geoms=rw_shapes, crs=proj,
                  edgecolor='blue', facecolor='blue')
ax.set_extent((-125, -118, 32, 49.2))
plt.show()
