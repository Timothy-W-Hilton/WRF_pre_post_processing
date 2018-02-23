"""functions useful for interpreting WRF grid coordinates
"""

import numpy as np
import netCDF4
from wrf import getvar
from shapely.geometry import Polygon
import geopandas as gp

def get_wrf_latlon(wrf_file, lonvar='XLONG', latvar='XLAT'):
    """read latitude and longitude from WRF input or output

    ARGS:
    wrf_file (str): full path to the WRF file
    lonvar (str): longitude variable in the netCDF file
    latvar (str): latitude variable in the netCDF file

    RETURNS:
    2-tuple of 2-D arrays (lon, lat)
    """
    nc = netCDF4.Dataset(wrf_file)
    lat = getvar(nc, latvar)
    lon = getvar(nc, lonvar)
    nc.close()
    return(lon, lat)


def get_cell_corner_coords(fname_wrf):
    """calculate 4 corner coordinates for every WRF grid cell

    Calculate the grid cell corners from in information in a WRF input
    or output file. fname_wrf must contain the netCDF variables
    XLONG_U, XLONG_V, XLAT_U, XLAT_V

    ARGS:
    fname_wrf (str): full path to a WRF data file.

    RETURNS:
    tuple of arrays: (lon_LL, lat_LL, lon_UL, lat_UL, lon_UR, lat_UR,
                      lon_LR, lon_LR)

    """
    lon_u, lat_u = get_wrf_latlon(fname_wrf,
                                  'XLONG_U', 'XLAT_U')
    lon_v, lat_v = get_wrf_latlon(fname_wrf,
                                  'XLONG_V', 'XLAT_V')
    # Calculate lower left, upper left, upper right, lower right cell
    # corner coordinates. For WRF grid illustration see
    # https://www.ncl.ucar.edu/Applications/Images/wrf_debug_3_lg.png
    lat_LL = lat_v.data[:-1, :]
    lon_LL = lon_u.data[:, :-1]

    lat_UL = lat_v.data[1:, :]
    lon_UL = lon_u.data[:, :-1]

    lat_UR = lat_v.data[1:, :]
    lon_UR = lon_u.data[:, 1:]

    lat_LR = lat_v.data[:-1, :]
    lon_LR = lon_u.data[:, 1:]

    return((lon_LL, lat_LL, lon_UL, lat_UL, lon_UR, lat_UR, lon_LR, lat_LR))


def WRF_cells_to_shapes_list(fname_wrf):
    """create a geopandas.geodataframe object for WRF corners
    """
    (lon_LL, lat_LL,
     lon_UL, lat_UL,
     lon_UR, lat_UR,
     lon_LR, lat_LR) = get_cell_corner_coords(fname_wrf)

    vertices_list = list(zip(
        list(zip(lon_LL.flatten(), lat_LL.flatten())),
        list(zip(lon_UL.flatten(), lat_UL.flatten())),
        list(zip(lon_UR.flatten(), lat_UR.flatten())),
        list(zip(lon_LR.flatten(), lat_LR.flatten()))))
    pgons_list = [Polygon(these_vertices) for these_vertices in vertices_list]
    idx = np.unravel_index(np.arange(lon_LL.size), lon_LL.shape)
    return(pgons_list, idx)

def pgons_2_gdf(pgons_list, idx):
    """create a geodataframe from a list of shapely polygons
    """
    gdf = gp.GeoDataFrame(geometry=pgons_list)
    gdf['x'] = idx[0]
    gdf['y'] = idx[1]
    return(gdf)
