"""functions useful for interpreting WRF grid coordinates
"""

import numpy as np
import netCDF4
from wrf import getvar
from shapely.geometry import Polygon
import geopandas as gp
from pyproj import Proj, transform

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


def _WRF_cells_to_shapes_list(fname_wrf, proj=None):
    """create a list of polygons representing WRF grid cells

    The WRF cells are represented as `shapely.geometry.Polygon` objects.

    ARGS:
    fname_wrf (str): full path to a WRF data file.
    proj (`proj4.Proj`): optional projection to apply to cell corner
       coordinates

    RETURNS
    list of `shapely.geometry.Polygon` objects reprensting WRF grid cells
    """
    (lon_LL, lat_LL,
     lon_UL, lat_UL,
     lon_UR, lat_UR,
     lon_LR, lat_LR) = get_cell_corner_coords(fname_wrf)

    if proj is not None:
        latlonProj = Proj(init='epsg:4326')
        x_LL, y_LL = transform(latlonProj, proj,
                               lon_LL, lat_LL)
        x_UL, y_UL = transform(latlonProj, proj,
                               lon_UL, lat_UL)
        x_UR, y_UR = transform(latlonProj, proj,
                               lon_UR, lat_UR)
        x_LR, y_LR = transform(latlonProj, proj,
                               lon_LR, lat_LR)
    else:
        x_LL, y_LL = (lon_LL, lat_LL)
        x_UL, y_UL = (lon_UL, lat_UL)
        x_UR, y_UR = (lon_UR, lat_UR)
        x_LR, y_LR = (lon_LR, lat_LR)

    vertices_list = list(zip(
        list(zip(x_LL.flatten(), y_LL.flatten())),
        list(zip(x_UL.flatten(), y_UL.flatten())),
        list(zip(x_UR.flatten(), y_UR.flatten())),
        list(zip(x_LR.flatten(), y_LR.flatten()))))
    pgons_list = [Polygon(these_vertices) for these_vertices in vertices_list]
    idx = np.unravel_index(np.arange(lon_LL.size), lon_LL.shape)
    return(pgons_list, idx)

def _pgons_2_gdf(pgons_list, idx):
    """create a geodataframe from a list of shapely polygons and indices

    The indices index the list items to the WRF grid [x, y] coordinates.

    ARGS:
    pgons_list (list): list of `shapely.geometry.Polygon` objects
       reprensting WRF grid cells, as from `WRF_cells_to_shapes_list()`.
    idx (two-tuple): tuple of lists containing x and y indices of the
       WRF cell polygons in pgons_list

    RETURNS
    `geopandas.GeoDataFrame` object containing the WRF cell polygons
       in its geometry and the indices in columns 'x' and 'y'
    """
    gdf = gp.GeoDataFrame(geometry=pgons_list)
    gdf['x'] = idx[0]
    gdf['y'] = idx[1]
    return(gdf)

def WRF_cells_to_gdf(fname_wrf, proj=None):
    """make a geopandas dataframe of polygons representing WRF grid cells

    this is a wrapper for WRF_cells_to_shapes_list, pgons_2_gdf

    ARGS:
    fname_wrf (str): full path to a WRF data file.
    proj (`proj4.Proj`): optional projection to apply to cell corner
       coordinates

    RETURNS
    `geopandas.GeoDataFrame` object containing the WRF cell polygons
       in its geometry and the indices in columns 'x' and 'y'
    """
    pgons_list, idx = _WRF_cells_to_shapes_list(fname_wrf, proj)
    gdf = _pgons_2_gdf(pgons_list, idx)
    return(gdf)
