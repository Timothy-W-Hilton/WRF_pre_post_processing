import os
import xarray as xr

def yatir_to_xarray(fname, varname, groupname=None, timerange=None):
    """read Yatir forest plot_diff.var_diff saved netCDF file to xarray

    returns an xarray suitable for plotting with
    [GeoViews](http://geoviews.org/user_guide/)
    """
    ds_all = xr.open_dataset(fname)
    if groupname is not None:
        ds_grp = xr.open_dataset(fname, group=groupname)

    var = xr.Dataset(data_vars={varname:
                                (['time', 'x', 'y'],
                                 ds_grp.variables[varname][slice(timerange), ...])},
                     coords={'time': (['time'], ds_all.time[slice(timerange), ...]),
                             'lat': (['x', 'y'], ds_all.lon),
                             'lon': (['x', 'y'], ds_all.lat)})
    return(var)


if __name__ == '__main__':

    LHctl = yatir_to_xarray(os.path.join('/', 'Users', 'tim', 'work',
                                         'Data', 'SummenWRF', 'yatir',
                                         'LH_d03_yatirZ50.nc'),
                            varname='LH',
                            groupname='ctl',
                            timerange=10)
    LHytr = yatir_to_xarray(os.path.join('/', 'Users', 'tim', 'work',
                                         'Data', 'SummenWRF', 'yatir',
                                         'LH_d03_yatirZ50.nc'),
                            varname='LH',
                            groupname='yatirZ050',
                            timerange=10)


# LH = xr.open_dataset('/Users/tim/work/Data/SummenWRF/yatir/LH_d03_yatirZ50.nc')
# LHctl = xr.open_dataset('/Users/tim/work/Data/SummenWRF/yatir/LH_d03_yatirZ50.nc',
#                         group='ctl')
# ntime = 10
# LH_byhand = xr.Dataset(data_vars={'LH': (['time', 'x', 'y'], LHctl.LH[:ntime, ...])},
#                       coords={'time': (['time'], np.arange(ntime)), # np.arange(LHctl.time.size)),
#                              'lat': (['x', 'y'], LH.lat),
#                               'lon': (['x', 'y'], LH.lon)})
