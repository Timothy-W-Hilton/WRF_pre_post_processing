"""demonstrate redwood location -- WRF grid cross referencing
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os.path
import wrf_grid
import redwoods_shapes
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import wrf_redwoods_crossref

def draw_map(gdf_wrf_all, gdf_wrf_rw, gdf_rw):
    """map WRF grid & redwood locations; highlight WRF cells containing redwoods

    ARGS:
    gdf_wrf_all (GeoDataFrame): GeoDataFrame containing boundaries of
       each WRF grid cell as a geosequence of polygons
    gdf_wrf_rw (GeoDataFrame): GeoDataFrame containing the subset of
       gdf_wrf_all that contain redwoods
    gdf_rw (GeoDataFrame): GeoDataFrame containing Redwood
       locations as geosequence of polygons

    """
    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')
    ax.add_feature(states_provinces, edgecolor='black')
    ax.add_geometries(geoms=gdf_wrf_rw.geometry, crs=ccrs.Mercator.GOOGLE,
                      facecolor='#1b9e77', edgecolor='black', alpha=0.5,
                      label='Redwood range')
    ax.add_geometries(geoms=gdf_wrf_all.geometry, crs=ccrs.Mercator.GOOGLE,
                      facecolor='white', edgecolor='black', alpha=0.3,

                      label='WRF cells with redwoods')
    ax.add_geometries(geoms=gdf_rw.geometry, crs=ccrs.Mercator.GOOGLE,
                      facecolor='#d95f02', edgecolor='black', alpha=0.5,
                      label='WRF cells with redwoods')
    bounds = gdf_wrf.total_bounds
    ax.set_extent(bounds[(0, 2, 1, 3), ], crs=ccrs.Mercator.GOOGLE)

    # make proxy artists for a legend
    redwoods = mpatches.Rectangle((0, 0), 1, 1, facecolor='#d95f02')
    wrf_cells = mpatches.Rectangle((0, 0), 1, 1, facecolor='white',
                                      edgecolor='black', alpha=0.3)
    wrf_redwoods = mpatches.Rectangle((0, 0), 1, 1, facecolor='#1b9e77',
                                   edgecolor='black', alpha=0.5)
    labels = ['Redwoods range',
              'WRF grid',
              'WRF cells with redwoods']
    plt.legend([redwoods, wrf_cells, wrf_redwoods], labels,
               loc='lower left', bbox_to_anchor=(0.025, -0.1), fancybox=True)

    plt.show()


if __name__ == "__main__":
    fname_wrf = os.path.join('/', 'Users', 'tim', 'work', 'Data',
                             'WRF_Driver',
                             'met_em.d01.2009-06-01_00:00:00.nc')
    fname_redwoods = os.path.join(
        '/', 'Users', 'tim', 'work', 'Data', 'Redwoods',
        'Redwood_SequoiaSempervirens_extentNorthAmerica',
        'data', 'commondata', 'data0', 'sequsemp.shp')

    (gdf_rw, gdf_wrf,
     wrf_cells_with_redwoods,
     redwoods_mask) = wrf_redwoods_crossref.main_wrapper(fname_redwoods, fname_wrf)

    draw_map(gdf_wrf, wrf_cells_with_redwoods, gdf_rw)
