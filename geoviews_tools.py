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

def build_geoviews_comparison(ds1, ds2):
    map_d02 = hv.Overlay((gf.land.options(scale='50m'),
                          gf.coastline.options(scale='50m'),
                          gf.borders.options(scale='50m'),
                          gv.Dataset(LHytr_d02).to(gv.QuadMesh,
                                                   groupby='time').opts(
                                                       xlim=(lat_d02.min() - pad02,
                                                             lat_d02.max() + pad02),
                                                       ylim=(lon_d02.min() - pad02,
                                                             lon_d02.max() + pad02)),
                          hv.Bounds((lat_d03.min(),
                                     lon_d03.min(),
                                     lat_d03.max(),
                                     lon_d03.max())).opts(color='blue')),
                         group=groupname,
                         label='d02')
    map_d03 = hv.Overlay((gf.land.options(scale='50m'),
                          gf.ocean.options(scale='50m'),
                          gf.coastline.options(scale='50m'),
                          gf.borders.options(scale='50m'),
                          gv.Dataset(LHytr_d03,
                                     group=groupname,
                                     label='d03').to(gv.QuadMesh,
                                                     groupby='time').opts(
                                                         xlim=(lat_d03.min() - pad03,
                                                               lat_d03.max() + pad03),
                                                         ylim=(lon_d03.min() - pad03,
                                                               lon_d03.max() + pad03)),
                          hv.Bounds((lat_d03.min(),
                                     lon_d03.min(),
                                     lat_d03.max(),
                                     lon_d03.max())).opts(color='blue')),
                         group=groupname,
                         label='d03')
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
