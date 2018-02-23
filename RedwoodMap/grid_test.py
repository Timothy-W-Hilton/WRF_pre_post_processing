import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os.path
import wrf_grid
import redwoods_shapes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from pyproj import Proj

def draw_map(gdf_wrf, gdf_rw):
    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')
    ax.add_feature(states_provinces, edgecolor='grey')
    ax.add_geometries(geoms=gdf_wrf.geometry, crs=ccrs.Mercator.GOOGLE,
                      facecolor=None, edgecolor='black', alpha=0.5,
                      label='Redwood range')
    ax.add_geometries(geoms=gdf_rw.geometry, crs=ccrs.Mercator.GOOGLE,
                      facecolor='#d95f02', edgecolor='black', alpha=0.5,
                      label='WRF cells with redwoods')
    bounds = gdf_wrf.total_bounds
    ax.set_extent(bounds[(0, 2, 1, 3), ], crs=ccrs.Mercator.GOOGLE)

    # make two proxy artists for a legend
    redwoods = mpatches.Rectangle((0, 0), 1, 1, facecolor='#d95f02')
    wrf_cells = mpatches.Rectangle((0, 0), 1, 1, facecolor=None, alpha=0.5)
    labels = ['Redwoods range',
              'WRF cells with redwoods']
    plt.legend([redwoods, wrf_cells], labels,
               loc='lower left', bbox_to_anchor=(0.025, -0.1), fancybox=True)

    plt.show()


fname_wrf = os.path.join('/', 'Users', 'tim', 'work', 'Data', 'WRF_Driver',
                         'met_em.d02.2009-06-01_00:00:00.nc')
fname_redwoods = os.path.join('/', 'Users', 'tim', 'work', 'Data',
                              'Redwoods',
                              'Redwood_SequoiaSempervirens_extentNorthAmerica',
                              'data', 'commondata', 'data0', 'sequsemp.shp')

rwProj = Proj(init='epsg:3857')
gdf_wrf = wrf_grid.WRF_cells_to_gdf(fname_wrf, proj=rwProj)
gdf_rw = redwoods_shapes.redwood_range_to_gdf(fname_redwoods)


gdf_wrf['has_redwoods'] = False
for i in range(gdf_rw.size):
    gdf_wrf['has_redwoods'] = (
        gdf_wrf['has_redwoods'] |
        gdf_wrf.intersects(gdf_rw.iloc[i].geometry))
wrf_cells_with_redwoods = gdf_wrf[gdf_wrf['has_redwoods']]

# this determines the shape of the intersection, it takes a lot longer
# than the above loop, which merely determines whether or not each WRF
# grid cell overlaps with any redwoods polygon
# cells_that_have_redwoods = gpd.overlay(gdf_wrf, gdf_rw, how='intersection')
