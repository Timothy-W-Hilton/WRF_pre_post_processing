import os
import xarray as xr
import holoviews as hv
import geoviews as gv
import geoviews.feature as gf


def yatir_to_xarray(fname, varname, groupname=None, timerange=None):
    """read Yatir forest plot_diff.var_diff saved netCDF file to xarray

    returns an xarray suitable for plotting with
    [GeoViews](http://geoviews.org/user_guide/)
    """
    ds_all = xr.open_dataset(fname)
    if groupname is not None:
        ds_grp = xr.open_dataset(fname, group=groupname)

    var = xr.Dataset(data_vars={varname: (['time', 'x', 'y'],
                                          ds_grp.variables[varname][
                                              slice(timerange), ...])},
                     coords={'time': (['time'],
                                      ds_all.time[slice(timerange), ...]),
                             'lat': (['x', 'y'], ds_all.lon),
                             'lon': (['x', 'y'], ds_all.lat)},
                     attrs={'varname': varname,
                            'groupname': groupname})
    return(var)


def build_geoviews_comparison(ds1, ds2, pad=[0.1, 1.0]):
    lon1 = ds1.lon.values
    lat1 = ds1.lat.values
    lon2 = ds2.lon.values
    lat2 = ds2.lat.values
    map_d02 = hv.Overlay((gf.land.options(scale='50m'),
                          gf.coastline.options(scale='50m'),
                          gf.borders.options(scale='50m'),
                          gv.Dataset(ds1).to(gv.QuadMesh,
                                             groupby='time').opts(
                                                 xlim=(lat1.min() - pad[0],
                                                       lat1.max() + pad[0]),
                                                 ylim=(lon1.min() - pad[0],
                                                       lon1.max() + pad[0])),
                          hv.Bounds((lat2.min(),
                                     lon2.min(),
                                     lat2.max(),
                                     lon2.max())).opts(color='blue')),
                         group=ds1.groupname,
                         label=ds1.varname)
    map_d03 = hv.Overlay((gf.land.options(scale='50m'),
                          gf.ocean.options(scale='50m'),
                          gf.coastline.options(scale='50m'),
                          gf.borders.options(scale='50m'),
                          gv.Dataset(ds2,
                                     group=ds2.groupname,
                                     label='d03').to(gv.QuadMesh,
                                                     groupby='time').opts(
                                                         xlim=(lat2.min() - pad[1],
                                                               lat2.max() + pad[1]),
                                                         ylim=(lon2.min() - pad[1],
                                                               lon2.max() + pad[1])),
                          hv.Bounds((lat2.min(),
                                     lon2.min(),
                                     lat2.max(),
                                     lon2.max())).opts(color='blue')),
                         group=ds2.groupname,
                         label=ds2.varname)
    return(map_d02, map_d03)


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
