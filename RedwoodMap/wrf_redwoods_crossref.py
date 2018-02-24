"""determine which WRF cells contain redwoods
"""

import numpy as np
from . import redwoods_shapes
import geopandas as gpd
from pyproj import Proj

import os.path
from . import wrf_grid

def get_WRF_geodataframe(fname_wrf, proj=None):
    gdf_wrf = wrf_grid.WRF_cells_to_gdf(fname_wrf, proj=proj)
    return(gdf_wrf)

def get_RW_geodataframe(fname_redwoods):
    gdf_rw = redwoods_shapes.redwood_range_to_gdf(fname_redwoods)
    return(gdf_rw)

def get_wrf_cells_with_redwoods(gdf_wrf, gdf_rw):

    gdf_wrf['has_redwoods'] = False
    for i in range(gdf_rw.size):
        gdf_wrf['has_redwoods'] = (
            gdf_wrf['has_redwoods'] |
            gdf_wrf.intersects(gdf_rw.iloc[i].geometry))
    wrf_cells_with_redwoods = gdf_wrf[gdf_wrf['has_redwoods']]
    # this should determine the shape of the intersection, it takes a
    # lot longer than the above loop, which merely determines whether
    # or not each WRF grid cell overlaps with any redwoods polygon
    # UNTESTED
    # cells_that_have_redwoods = gpd.overlay(gdf_wrf, gdf_rw, how='intersection')
    return(gdf_wrf, wrf_cells_with_redwoods)

def get_redwood_mask(gdf_wrf, wrf_cells_with_redwoods):
    # add 1 to max indices to account for zero-based indices
    redwood_mask = np.zeros((gdf_wrf['x'].max() + 1,
                             gdf_wrf['y'].max() + 1),
                            dtype=bool)
    redwood_mask[wrf_cells_with_redwoods['x'],
                 wrf_cells_with_redwoods['y']] = True
    return(redwood_mask)

def main_wrapper(fname_redwoods, fname_wrf):
    rwProj = Proj(init='epsg:3857')
    gdf_rw = get_RW_geodataframe(fname_redwoods)
    gdf_wrf = get_WRF_geodataframe(fname_wrf, proj=rwProj)
    gdf_wrf, wrf_cells_with_redwoods = get_wrf_cells_with_redwoods(gdf_wrf,
                                                                   gdf_rw)
    redwoods_mask = get_redwood_mask(gdf_wrf, wrf_cells_with_redwoods)
    return(gdf_rw, gdf_wrf, wrf_cells_with_redwoods, redwoods_mask)

if __name__ == "__main__":
    fname_wrf = os.path.join('/', 'Users', 'tim', 'work', 'Data',
                             'WRF_Driver',
                             'met_em.d02.2009-06-01_00:00:00.nc')
    fname_redwoods = os.path.join(
        '/', 'Users', 'tim', 'work', 'Data', 'Redwoods',
        'Redwood_SequoiaSempervirens_extentNorthAmerica',
        'data', 'commondata', 'data0', 'sequsemp.shp')

    (gdf_rw, gdf_wrf,
     wrf_cells_with_redwoods,
     redwoods_mask) = main_wrapper(fname_redwoods, fname_wrf)
