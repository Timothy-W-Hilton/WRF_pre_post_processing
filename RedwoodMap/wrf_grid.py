"""functions useful for interpreting WRF grid coordinates
"""

import netCDF4
from wrf import getvar

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
    lat_LL = lat_v[:-1, :]
    lon_LL = lon_u[:, :-1]

    lat_UL = lat_v[1:, :]
    lon_UL = lon_u[:, :-1]

    lat_UR = lat_v[1:, :]
    lon_UR = lon_u[:, 1:]

    lat_LR = lat_v[:-1, :]
    lon_LR = lon_u[:, 1:]

    return((lon_LL, lat_LL, lon_UL, lat_UL, lon_UR, lat_UR, lon_LR, lat_LR))
