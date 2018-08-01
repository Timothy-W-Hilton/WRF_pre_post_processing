import cartopy.feature as cfeature
import numpy as np
import os.path
import prism_tools
import pickle

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs


parse_data = False


def plot_interpolated(pts):
    """
    """

    # set up axes
    fig = plt.figure(figsize=(11, 6))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.05)
    ax = [fig.add_subplot(gs[0, n], projection=ccrs.PlateCarree())
          for n in range(2)]
    for this_ax in ax:
        this_ax.set_extent([pts.lon.min() - 1.0, pts.lon.max() + 1.0,
                            pts.lat.min() - 1.0, pts.lat.max() + 1.0],
                           crs=ccrs.PlateCarree())
        this_ax.add_feature(cfeature.LAND)
        this_ax.add_feature(cfeature.OCEAN)
    # plot data
    long, latg = np.meshgrid(pts.lon, pts.lat)
    cs = ax[0].pcolormesh(long, np.flipud(latg), pts.data[0, ...])
    ax[0].set_title('original PRISM')
    cs = ax[1].pcolormesh(long, np.flipud(latg), pts.data[0, ...])
    ax[1].set_title('interpolated NN')
    # colorbar
    ax = fig.add_subplot(gs[0, 2])
    ax.set_title(pts.varname)
    plt.colorbar(cs, cax=ax)
    return(fig)

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
    # new_data = pts.interpolate(lon, lat, method='NN')
    plot_interpolated(pts)
