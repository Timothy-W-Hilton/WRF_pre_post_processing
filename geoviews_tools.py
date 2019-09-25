import numpy as np
import os
import xarray as xr
import holoviews as hv
from holoviews import opts
import geoviews as gv
import geoviews.feature as gf

from adjust_WRF_inputs import km_to_yatir


def yatir_mask(xarr, d_threshold=5):
    """add a variable identifying WRF cells in Yatir parameterization

    The new variable is boolean, True for cells in Yatir, False otherwise.

    ARGS:
    xarr [xarray.Dataset]: dataset containing Yatir WRF data.  Must
       contain variables 'lon' and 'lat'
    d_threshold [int]: distance threshold (in kilometers) to be
       considered "in Yatir".  Default is 5 km.
    """
    d_to_yatir = km_to_yatir(xarr['lat'].values, xarr['lon'].values)
    xarr['yatir_mask'] = (('x', 'y'), d_to_yatir < d_threshold)
    return(xarr)


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
    var = yatir_mask(var)
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


def xrvar_to_hvdim(var):
    """create a HoloViews dimension from an xarray variable

    Assumes that the xarray variable has attributes name, description,
    and units

    RETURNS:
        hv.Dimension object

    """
    try:
        return(hv.Dimension(var.name, label=var.description, unit=var.units))
    except AttributeError:
        print('var must have attributes name, description, and units')
        raise


def overlay_roughness_realization_timeseries(dsxr, varname):
    """create a HoloViews Overlay showing one variable's timeseries

    to create a timeseries from two-dimensional data, the data are
    averaged spatially before plotting

    RETURNS:
      hvOverlay object containing timeseries for all WRF runs in the experiment

    """
    df_var = dsxr.mean(dim=['west_east', 'south_north'],
                       skipna=True)[varname].to_dataframe()
    vert_dim = xrvar_to_hvdim(dsxr[varname])
    tbl_list = [hv.Table(df[1], label=df[0]) for
                df in df_var.groupby('WRFrun')]
    crv_list = [hv.Curve(tbl, kdims='Time', vdims=vert_dim) for
                tbl in tbl_list]
    OLopts = opts.Overlay(frame_width=690, legend_position='right')
    label = 'roughness length experiments: {}'.format(
        dsxr[varname].description)
    hvol = hv.Overlay(crv_list).opts(OLopts).relabel(label)
    return(hvol)


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
