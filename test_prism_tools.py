import cartopy.feature as cfeature
import numpy as np
import os.path
import prism_tools
import pickle

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs

import interpolator

parse_data = False


def update(event, pts, lon, lat, ax):
    if pts.tidx < pts.data.shape[0]:
        pts.tidx = pts.tidx + 1
    else:
        pts.tidx = 0
    long, latg = np.meshgrid(pts.lon, pts.lat)
    #clear frame
    # ax[0].clear()
    # ax[1].clear()
    ax[0].pcolormesh(long, latg, pts.data[pts.tidx, ...],
                     vmin=0, vmax=30)
    ax[1].pcolormesh(lon, lat, pts.data_interp[pts.tidx, ...],
                     vmin=0, vmax=30)
    plt.draw() #redraw


def plot_interpolated(pts, lon, lat):
    """
    """

    # set up axes
    fig = plt.figure(figsize=(11, 6))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.05)
    ax = [fig.add_subplot(gs[0, n], projection=ccrs.PlateCarree())
          for n in range(2)]
    for this_ax in ax:
        this_ax.set_extent([lon.min() - 1.0, lon.max() + 1.0,
                            lat.min() - 1.0, lat.max() + 1.0],
                           crs=ccrs.PlateCarree())
        this_ax.coastlines(resolution='50m', color='black', linewidth=1)
    # plot data
    long, latg = np.meshgrid(pts.lon, pts.lat)
    cs = ax[0].pcolormesh(long, latg, pts.data[pts.tidx, ...],
                          vmin=0, vmax=30)
    ax[0].set_title('original PRISM')
    cs = ax[1].pcolormesh(lon, lat, pts.data_interp[pts.tidx, ...],
                          vmin=0, vmax=30)
    ax[1].set_title('interpolated NN')
    # colorbar
    ax_cb = fig.add_subplot(gs[0, 2])
    ax_cb.set_title(pts.varname)
    plt.colorbar(cs, cax=ax_cb)
    return(fig, ax)


if __name__ == "__main__":
    prism_dir = os.path.join('/', 'Users',
                             'tim', 'work', 'Data', 'PRISM')
    fname_pickle = "pts.pickle"
    if parse_data:

        pts = prism_tools.PRISMTimeSeries(
            os.path.join(prism_dir,
                         'PRISM_tmean_stable_4kmD1_20090601_20090630_asc.zip'),
            varname='tmean', varunits='C')
        pts.from_zip_archive()
        pfile = open(fname_pickle, 'wb')
        pickle.dump(pts, pfile)
        pfile.close()
    else:
        pfile = open(fname_pickle, 'rb')
        pts = pickle.load(pfile)
        pfile.close()

    lon, lat = prism_tools.read_WRF_latlon(
        os.path.join(prism_dir, 'WRF_d02_latlon.nc'))
    pts.interpolate(lon, lat, method='NN')
    pts.tidx = 0
    fig, ax = plot_interpolated(pts, lon, lat)
    fig.canvas.mpl_connect('button_press_event',
                           lambda event: update(1, pts, lon, lat, ax))
    plt.show()
