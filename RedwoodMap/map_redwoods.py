"""script to plot redwood location polygons to a map

modified from https://gis.stackexchange.com/questions/259831/shapefiles-wont-display-with-cartopy
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

fig = plt.figure(figsize=(5, 5))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution='10m')
countries = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_boundary_lines_land',
    scale='10m',
    facecolor='none')
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')
oceans = cfeature.NaturalEarthFeature(
    category='physical',
    name='Ocean',
    scale='50m',
    facecolor='#87ceeb')
ax.add_feature(countries, edgecolor='grey')
ax.add_feature(states_provinces, edgecolor='grey')
ax.add_feature(oceans, facecolor='#87ceeb')
ax.add_geometries(geoms=rw_shapes, crs=proj,
                  edgecolor='#d95f02', facecolor='#d95f02')
ax.set_extent((-124.5, -114, 31.5, 42.5))
long_marine_lab = (-122.064791, 36.949690)
ax.annotate('Long Marine Lab',
            xy = long_marine_lab,
            xytext = (-122.064791 - 2, 36.949690 - 2),
            xycoords = 'data',
            arrowprops={'arrowstyle': '-|>'},
            transform=ccrs.PlateCarree())
plt.show()
fig.savefig('/Users/tim/Desktop/test.pdf')
