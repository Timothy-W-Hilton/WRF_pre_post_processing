import numpy as np
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
                             'lon': (['x', 'y'], ds_all.lat),
                             'x': (['x'], np.arange(
                                 ds_grp.variables[varname].shape[1])),
                             'y': (['y'], np.arange(
                                 ds_grp.variables[varname].shape[2]))},
                     attrs={'varname': varname,
                            'groupname': groupname,
                            'units': ds_all.units})
    return(var)


def get_quadmesh(data, pad):
    qm = (gv.Dataset(data).to(gv.QuadMesh, groupby='time').opts(
        xlim=(data.lat.values.min() - pad,
              data.lat.values.max() + pad),
        ylim=(data.lon.values.min() - pad,
              data.lon.values.max() + pad)))
    return(qm)


def build_geoviews_comparison(ds_d02, ds_d03, pad=[1.0, 0.1]):
    """build four maps showing a WRF variable for [control, Yatir] and [d02, d03]
    """
    # overlay land, oceans, coasts, political borders
    # map_background = (gf.land.options(scale='50m'),
    #                   gf.ocean.options(scale='50m'),
    #                   gf.coastline.options(scale='50m'),
    #                   gf.borders.options(scale='50m'))
    map_background = (gf.coastline.options(scale='50m'),
                      gf.borders.options(scale='50m'))
    # make bounding box for inner-most nested domain
    d03_bb = hv.Bounds((ds_d03.lat.values.min(),
                        ds_d03.lon.values.min(),
                        ds_d03.lat.values.max(),
                        ds_d03.lon.values.max())).opts(color='blue')
    map_d02 = hv.Overlay((get_quadmesh(ds_d02, pad[0]), ) +
                         map_background +
                         (d03_bb, ),
                         group=ds_d02.groupname,
                         label=ds_d02.varname)
    map_d03 = hv.Overlay((get_quadmesh(ds_d03, pad[1]), ) +
                         map_background +
                         (d03_bb, ),
                         group=ds_d03.groupname,
                         label=ds_d03.varname)
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
