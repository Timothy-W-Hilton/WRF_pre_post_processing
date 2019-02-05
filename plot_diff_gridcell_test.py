import datetime
import os
from plot_diff import var_diff, VarDiffPlotter
from wrf import ll_to_xy, ll_points
import numpy as np
import netCDF4

DOMAIN = 2
SanFrancisco = (37.7707405, -122.450006)
SantaCruz = (36.974474, -122.028986)

wrffile = os.path.join('/', 'Users', 'tim', 'work', 'Data',
                       'SummenWRF',
                       'metem_2005_ctl_NCEPDOE_d02_2005-06-01_00:00:00.nc')


def get_wrf_xy(coords):
    """return XY coordinates for a tuple of lat, lon pairs

    calculates XY for Summen WRF domain 2
    """
    nc = netCDF4.Dataset(wrffile, 'r')
    lat, lon = zip(*coords)
    xy = ll_to_xy(nc, lat, lon, stagger='u')
    nc.close()
    return(xy)


def compare_coords():
    """compare wrf.ll_points() to XLONG, XLAT variables
    """
    nc = netCDF4.Dataset(wrffile, 'r')
    xlat = nc.variables['XLAT_U'][...].squeeze()
    xlong = nc.variables['XLONG_U'][...].squeeze()
    ll_points_list = ll_points(xlat, xlong)
    return(ll_points_list, xlat, xlong)


if __name__ == "__main__":

    varname = 'fogpresent'

    read_data = False
    vd = var_diff(
        ncfile=os.path.join(
            '/Users/tim/work/Data/SummenWRF/',
            '{varname}_d{DOMAIN:02d}_CLM_nourban.nc'.format(
                varname=varname, DOMAIN=DOMAIN)))
    # for k in vd.data.keys():
    #     vd.data[k] = vd.data[k] * 100.0
    if vd.p is None:
        vd.get_significance_mask(significance=0.95, adj_autocorr=True)
    pfx = 'gridtest'
    # for this_series in ['all_tstamps', 'time_avg']:

    t_end = 1
    pfx = pfx + '_timeavg'
    vd.aggregate_time(time_avg=True)
    time_title_str = 'June 2005'

    for k in vd.data.keys():
        vd.data[k][...] = np.ones(vd.data[k].shape)

    xy = get_wrf_xy((SanFrancisco, SantaCruz))
    vd.data['ctl'][0, xy.data[1, 0], xy.data[0, 0]] = 50
    vd.data['no_urban_CLM'][0, xy.data[1, 1], xy.data[0, 1]] = 100

    for this_t in range(0, t_end):  #
        plotter = VarDiffPlotter(vd, t_idx=this_t, layer=0,
                                 domain=DOMAIN,
                                 pfx=pfx,
                                 savedir='.',
                                 time_title_str=time_title_str)
        if DOMAIN == 1:
            cb_orientation = 'horizontal'
        else:
            cb_orientation = 'vertical'
        fig, ax = plotter.plot(cb_orientation=cb_orientation,
                           vmin=0,
                           vmax=100,
                           mask=None)  #vd.p > 0.05)
