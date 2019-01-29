"""plot WRF dominant land use on a map
"""

import numpy as np
import netCDF4

# from matplotlib.cm import get_cmap
# from matplotlib.figure import Figure
import cartopy.crs as ccrs
import datetime

from map_tools_twh.map_tools_twh import Fig
from map_tools_twh.map_tools_twh import get_IGBP_modMODIS_21Category_PFTs_cmap
from map_tools_twh.map_tools_twh import get_IGBP_modMODIS_21Category_PFTs_table
from timutils.colormap_nlevs import setup_colormap
from plot_diff import var_diff, VarDiffPlotter

bbox_SFBay = (-123.0, -121.0, 36.2, 39.0)
bbox_LA = (-120.0, -116.0, 32.0, 35.0)


class LU_vardiff(var_diff):
    """subclass of var_diff to deal with land use

    land use is categorical data, not numeric, so tabulating
    differencing is, well, different :)

    """


def plot_init(lon, lat):

    # initialize figure, axes
    nplots = 6
    fig = Fig(figsize=(15, 15))
    ax = [None] * nplots
    for axidx, axspec in enumerate(range(321, 321 + nplots)):
        if axidx < 10:
            prj = ccrs.PlateCarree()
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
                   lat.min(), lat.max()))
    ax.coastlines(resolution='10m', color='black')
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
    axes_right_side = (0, 2, 4)
    axes_left_side = (1, 3, 5)
    for this_ax in axes_right_side:
        cm = plot_landuse(ax[this_ax], lon, lat, luidx['ctl'])
    for this_ax in axes_left_side:
        cm = plot_landuse(ax[this_ax], lon, lat, luidx['deurb'])

    ax[0].set_title('control')
    ax[1].set_title('deurbanized')
    ax[2].set_extent(bbox_SFBay, crs=ccrs.PlateCarree())
    ax[3].set_extent(bbox_SFBay, crs=ccrs.PlateCarree())
    ax[4].set_extent(bbox_LA, crs=ccrs.PlateCarree())
    ax[5].set_extent(bbox_LA, crs=ccrs.PlateCarree())

    cbar = fig.colorbar(cm,
                        ax=ax,
                        cmap=cmap,
                        extend='neither',
                        ticks=np.linspace(0.5, 21.5, 21))
    cbar.set_ticks(np.linspace(0.5, 21.5, 21))
    cbar.set_ticklabels(list(ctable['long_name']))
    # for axidx in range(6):
    #     print('{} coastlines'.format(axidx))
    #     ax[axidx].coastlines('10m', color='black')

    fig.savefig(fname='/tmp/foo.png')
