"""plot WRF dominant land use on a map
"""

import numpy as np
import numpy.ma as ma
import datetime
import os
import netCDF4
import glob
from xarray import DataArray
from wrf import getvar, extract_times, to_np, ALL_TIMES

# from matplotlib.cm import get_cmap
# from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from map_tools_twh.map_tools_twh import CoastalSEES_WRF_prj
from map_tools_twh.map_tools_twh import Fig
from map_tools_twh.map_tools_twh import get_IGBP_modMODIS_21Category_PFTs_cmap
from map_tools_twh.map_tools_twh import get_IGBP_modMODIS_21Category_PFTs_table
from timutils.colormap_nlevs import setup_colormap
from plot_diff import var_diff, wrf_var


class LU_vardiff(var_diff):
    """subclass of var_diff to deal with land use

    land use is categorical data, not numeric, so tabulating
    differencing is, well, different :)

    """


def plot_init(lon, lat):

    # initialize figure, axes
    nplots = 4
    fig = Fig(figsize=(8, 8))
    ax = [None] * nplots
    for axidx, axspec in enumerate(range(221, 221 + nplots)):
        if axidx < 3:
            prj = CoastalSEES_WRF_prj()
        else:
            prj = None
        ax[axidx] = fig.add_subplot(axspec,
                                    projection=prj)
        # ax[axidx].set_extent((lon.min(), lon.max(),
        #                       lat.min(), lat.max()))
    return(fig, ax)


def plot_landuse(ax, lon, lat, data):
    cmap, norm = setup_colormap(
        0,
        21,
        nlevs=21,
        cmap=get_IGBP_modMODIS_21Category_PFTs_cmap())
    cm = ax.pcolormesh(lon,
                       lat,
                       data,
                       cmap=cmap,
                       norm=norm)
    ax.set_extent((lon.min(), lon.max(),
                   lat.min(), lat.max()),
                  crs=CoastalSEES_WRF_prj())
    return(cm)

if __name__ == "__main__":
    wrfin = {'ctl': '/Users/tim/work/Data/WRF_Driver/wrfinput_d02_ctl',
             'deurb': '/Users/tim/work/Data/WRF_Driver/wrfinput_d02_deurbanized'}
    lufrac = {}
    luidx = {}
    for k in wrfin.keys():
        nc = netCDF4.Dataset(wrfin[k], 'r')
        lufrac[k] = nc.variables['LANDUSEF'][...].squeeze()
        luidx[k] = nc.variables['LU_INDEX'][...].squeeze()
        lon = nc.variables['XLONG'][...].squeeze()
        lat = nc.variables['XLAT'][...].squeeze()
        nc.close()

    ctable = get_IGBP_modMODIS_21Category_PFTs_table()
    cmap = get_IGBP_modMODIS_21Category_PFTs_cmap()

    fig, ax = plot_init(lon, lat)
    cm = plot_landuse(ax[0], lon, lat, luidx['ctl'])
    cm = plot_landuse(ax[1], lon, lat, luidx['deurb'])
    cbar = fig.colorbar(cm,
                        ax=ax,
                        cmap=cmap,
                        extend='neither',
                        ticks=np.linspace(0.5, 21.5, 21))
    # tick_locs = (np.arange(cmap.N) + 0.5)*(cmap.N-1)/cmap.N
    # cbar.set_ticks(tick_locs)
    cbar.set_ticks(np.linspace(0.5, 21.5, 21))
    cbar.set_ticklabels(list(ctable['long_name']))
    fig.savefig(fname='/tmp/foo.png')
