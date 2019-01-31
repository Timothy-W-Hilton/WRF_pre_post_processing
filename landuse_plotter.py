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

    def __init__(self, fname_A=None, fname_B=None, label_A=None, label_B=None):
        super(LU_vardiff, self).__init__(
            fname_A=fname_A, fname_B=fname_B,
            label_A=label_A, label_B=label_B,
            varname='LANDUSEF')

    def read_files(self):
        for k in self.data.keys():
            nc = netCDF4.Dataset(self.fnames[k], 'r')
            self.data[k] = nc.variables['LANDUSEF'][...].squeeze()
            self.lon = nc.variables['XLONG'][...].squeeze()
            self.lat = nc.variables['XLAT'][...].squeeze()
            nc.close()
        self.time = (datetime.datetime(2005, 6, 1, 0, 0, 0), )

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
        lon = nc.variables['XLONG_U'][:, :, :-1].squeeze()
        lat = nc.variables['XLAT_V'][:, :-1, :].squeeze()
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

    fig.savefig(fname='/tmp/land_use_dominant.png')

    vd_LUfrac = LU_vardiff(fname_A=wrfin['ctl'], fname_B=wrfin['deurb'],
                           label_A='ctl', label_B='deurb')
    for this_pft in range(len(ctable)):
        vd_LUfrac.read_files()
        for k in vd_LUfrac.data.keys():
            vd_LUfrac.data[k] = np.expand_dims(
                vd_LUfrac.data[k][this_pft, ...],
                0)
        plotter = VarDiffPlotter(
            vd_LUfrac,
            t_idx=0,
            layer=0,
            domain=2,
            pfx='{:02d}-{}'.format(
                ctable['PFTnum'][this_pft],
                ctable['long_name'][this_pft].replace('/', '')),
            savedir='.',
            time_title_str='foo')
        fig = plotter.plot(vmin=0.0, vmax=1.0)
