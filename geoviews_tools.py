import numpy as np
import os
import socket
import xarray as xr
import holoviews as hv
from holoviews import opts
import geoviews as gv
import geoviews.feature as gf
import pandas as pd

from map_tools_twh.map_tools_twh import get_IGBP_modMODIS_21Category_PFTs_table
from adjust_WRF_inputs import km_to_yatir


def daily_cycle_mean_overlay_layout(df_obs, ds_WRF, varname, dims=None):
    """plot daily cycles for Yatir obs and WRF in two panels

    Left panel shows WRF desert cells, right panel shows WRF Yatir cells.

    ARGS:
    ds_dc (xarray.Dataset): xarray dataset ("ds") containing mean
       daily cycles ("dc".  Must have dimensions [area, hour, WRFrun].

    RETURNS:
    holoviews.HoloMap containing the plots
    """

    if dims is None:
        hour_dim = 'hour'
        var_dim = varname
    else:
        hour_dim = dims['hour']
        var_dim = dims[varname]

    # place curves for the requested variable in a holomap indexed by
    # area (desert, Yatir), "WRFrun" (control, Yatir, observations)
    ds_dc = combine_yatir_obs_WRF(df_obs, ds_WRF, varname)
    hm_cv = hv.Dataset(ds_dc).to(hv.Curve,
                                 kdims=[hour_dim],
                                 vdims=[var_dim])
    # hm_cv = hv.Curve(hv.Dataset(ds_dc),
    #                  kdims=hour_dim,
    #                  vdims=var_dim)
    opts_curve = opts.Curve(width=400)
    opts_ovly = opts.Overlay(legend_position='right')
    # place the holomap of curves into a holomap containing an
    # NdOverlay, indexed by WRFrun
    hm_ndovly = hm_cv.opts(opts_curve).overlay('WRFrun')
    # rearrange the holomap into an NdLayout
    ndovly = hm_ndovly.opts(opts_ovly).layout('area')

    hm_ndovly = hm_cv.overlay('WRFrun')
    # rearrange the holomap into an NdLayout
    ndovly = hm_ndovly.layout('area')

    return(ndovly)


def combine_yatir_obs_WRF(df_obs, ds_WRF, varname):
    """combine WRF output and Yatir observations into an xarray dataset
    """
    # create xarray DataArray from pandas dataframe column for varname and time
    da_obs = xr.DataArray(
        data=df_obs[varname].values[np.newaxis, np.newaxis, :],
        coords={'hour': (['hour'], df_obs.reset_index()['hour']),
                'area': (['area'], ['yatir']),
                'WRFrun': (['WRFrun'], ['Yatir obs'])},
        dims=('area', 'WRFrun', 'hour')
    )
    WRF_with_obs = xr.concat((da_obs, ds_WRF[varname]), dim='WRFrun')
    # not sure why I need this, but seems to be necessary to give the
    # DataArray a name
    WRF_with_obs = WRF_with_obs.rename(varname)
    return(WRF_with_obs)

def merge_yatir_fluxes_landuse():
    """merge WRF fluxes and landuse into single xarray dataset
    """
    cscratch_path = os.path.join('/', 'global', 'cscratch1', 'sd',
                                 'twhilton', 'yatir_output_collected')
    ctlday = WRF_daily_daylight_avg(os.path.join(cscratch_path,
                                                 'ctl_run_d03_diag.nc'))
    ytrday = WRF_daily_daylight_avg(os.path.join(cscratch_path,
                                                 'yatir_run_d03_diag.nc'))
    landuse_data = yatir_landuse_to_xarray()

    ytrday = ytrday.assign(
        {'LU_INDEX':
         landuse_data['d03'].sel(WRFrun='ytr')['LU_INDEX']})
    ctlday, ytrday = (this_dataset.assign(
        {'LU_INDEX':
         landuse_data['d03'].sel(WRFrun=this_key)['LU_INDEX'],
         'LANDUSEF':
         landuse_data['d03'].sel(WRFrun=this_key)['LANDUSEF']})
                      for (this_dataset, this_key) in
                      zip((ctlday, ytrday), ('ctl', 'ytr')))
    ds_diff =  (ctlday - ytrday).assign_coords({'WRFrun': 'control - Yatir'})
    return(ctlday, ytrday, ds_diff)


def yatir_mask(xarr, d_threshold=5,
               lonvar='lon', latvar='lat',
               hdim1='x', hdim2='y'):
    """add a variable identifying WRF cells in Yatir parameterization

    The new variable is boolean, True for cells in Yatir, False otherwise.

    ARGS:
    xarr [xarray.Dataset]: dataset containing Yatir WRF data.  Must
       contain variables 'lon' and 'lat'
    d_threshold [int]: distance threshold (in kilometers) to be
       considered "in Yatir".  Default is 5 km.
    lonvar (str): name of the longitude variable in xarr
    latvar (str): name of the latitude variable in xarr
    hdim1 (str): name of the first horizontal dimension in xarr
    hdim2 (str): name of the second horizontal dimension in xarr
    """

    d_to_yatir = km_to_yatir(xarr[latvar].values, xarr[lonvar].values)
    xarr['yatir_mask'] = ((hdim1, hdim2), d_to_yatir < d_threshold)
    return(xarr)


def get_data_file(data_paths):

    data = {k: xr.open_dataset(v).set_index(Time='XTIME')
            for k, v in data_paths.items()}

    datasets = list(data.values())

    dsxr = xr.concat(datasets,
                     dim=pd.Index(list(data.keys()), name='WRFrun'))
    # remove Time variation from latitude, longtitude arrays.  They
    # are in fact constant in time, and keeping the non-varying time
    # dimension around makes plotting more painful.
    for this_var in ['XLAT', 'XLONG',
                     'XLAT_U', 'XLONG_U',
                     'XLAT_V', 'XLONG_V']:
        dsxr[this_var] = dsxr[this_var].sel(Time=dsxr['Time'][1])
        dsxr[this_var] = dsxr[this_var].sel(WRFrun=dsxr['WRFrun'][0])
    # put in coordinate values for north/south and east/west
    # directions.  This seems to be necessary for GeoViews to properly
    # interpret the dataset when cast as a GeoViews Dataset.
    dsxr = dsxr.assign_coords(south_north=range(dsxr.south_north.size))
    dsxr = dsxr.assign_coords(west_east=range(dsxr.west_east.size))
    # get a mask for which pixels are parameterized to Yatir Forest
    dsxr = yatir_mask(dsxr,
                      latvar='XLONG', lonvar='XLAT',
                      hdim1='west_east', hdim2='south_north')
    # calculate differences between WRF runs
    # dsxr['d'] = dsxr.diff(dim='WRFrun').sel(WRFrun=
    # cast the xarray dataset to a GeoViews dataset
    dsgv = gv.Dataset(data=dsxr,
                      kdims=['Time', 'WRFrun', 'west_east', 'south_north'],
                      vdims=['XLAT', 'XLONG', 'HFX', 'LH', 'U', 'V', 'W'])
    return(dsgv, dsxr)


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


def yatir_WRF_to_xarray(fname):
    """read Yatir forest WRF output to xarray

    returns an xarray suitable for plotting with
    [GeoViews](http://geoviews.org/user_guide/)
    """
    ds = xr.open_dataset(fname)
    for dim in ['XLAT', 'XLONG']:
        ds[dim] = ds[dim].sel(Time=0)
    return(ds)


def WRF_yatir_desert_timeseries(fname):
    """calculate mean time series for Yatir, desert WRF cells

    TODO: implement )
    """


def WRF_daily_daylight_avg(fname):
    """read Yatir forest WRF output to xarray containing daily means

    returns an xarray suitable for plotting with
    [GeoViews](http://geoviews.org/user_guide/)
    """
    ds = yatir_WRF_to_xarray(fname)
    is_daytime = ds['SWDOWN'] > 0.1
    ds_day_mean = ds.where(is_daytime, drop=True).groupby('XTIME.hour').mean(keep_attrs=True)
    return(ds_day_mean)


def define_dims(ds):
    dims = {k: hv.Dimension(k, label=v.attrs['description'], unit=v.attrs['units'])
        for k, v in ds.data_vars.items()}
    dims['date'] = hv.Dimension('XTIME', label='date', unit='UTC')
    dims['hour'] = hv.Dimension('hour', label='hour of day', unit='UTC')
    dims['lon'] = hv.Dimension('XLONG', label='longitude', unit='deg E')
    dims['lat'] = hv.Dimension('XLAT', label='latitude', unit='deg N')
    return(dims)

# def yatir_WRF_diff(ds1, ds2, varname):
#     """calculate difference in a variable between ds1 and ds2
#     """
#     ds_diff = xr.concat(

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


def overlay_roughness_realization_timeseries(dsxr, varname, mask=None):
    """create a HoloViews Overlay showing one variable's timeseries

    to create a timeseries from two-dimensional data, the data are
    averaged spatially before plotting

    RETURNS:
      hvOverlay object containing timeseries for all WRF runs in the experiment

    """
    if mask is not None:
        dsxr = dsxr.where(mask)
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


def parse_yatir_EC_observations():
    """parse Yatir forest eddy covariance observations to a pandas data frame
    """
    if 'MacBook' in socket.gethostname():
        dir_path = os.path.join('/', 'Users', 'tim', 'work', 'Data',
                                'Yatir_Forest_Data')
    elif 'cori' in socket.gethostname():
        dir_path = os.path.join('/', 'project', 'projectdirs',
                                'm2319', 'Data', 'Yatir_Forest_Data')
    df_ytr = pd.read_csv(os.path.join(dir_path,
                                      'EFDC_L2_Flx_ILYat_2015_v03_30m.txt'),
                         na_values=[-9999])
    df_ytr['time'] = pd.to_datetime(df_ytr['TIMESTAMP_START'],
                                    format='%Y%m%d%H%M')
    df_ytr['hour'] = pd.DatetimeIndex(df_ytr['time']).hour
    return(df_ytr)



def yatir_landuse_to_xarray():
    """parse land use data for Yatir, Control runs for domains d02 and d03

    RETURNS:
      dict keyed by ['d02', 'd03']; values are xarray.DataSet objects
      containing land use data concatenated on new dimension WRFrun
    """
    if 'MacBook' in socket.gethostname():
        dir_path = os.path.join('/', 'Users', 'tim', 'work',
                                'Data', 'SummenWRF', 'yatir')
    elif 'cori' in socket.gethostname():
        dir_path = os.path.join('/', 'global', 'cscratch1', 'sd',
                                'twhilton', 'yatir_land_use')
    ctable = get_IGBP_modMODIS_21Category_PFTs_table()
    land_cat_names = list(ctable['long_name'])
    land_cat_names = [x if x != 'BareGroundTundra' else 'Yatir'
                      for x in land_cat_names]
    dict_runs = {}
    for WRFdomain in ['d02', 'd03']:
        dict_runs[WRFdomain] = xr.concat((xr.open_dataset(
            os.path.join(dir_path, 'land_data_{run}_{dom}.nc').format(
                run=WRFrun, dom=WRFdomain)).squeeze()
                                          for WRFrun in ['ctl', 'ytr']),
                                         dim='WRFrun')
        # remove the hyphen from z dimension name.  '-' also being an
        # operator messes up assign_coords()
        # dict_runs[WRFdomain] = dict_runs[WRFdomain].rename(
        #     {'z-dimension0021': 'zdimension0021'})

        # assign integral coordinate values to spatial coordinate
        # variables
        new_coords = {this_var: range(dict_runs[WRFdomain][this_var].size) for
                      this_var in ['south_north', 'west_east']}
        new_coords = {**new_coords,
                      'WRFrun': ['ctl', 'ytr'],
                      'z-dimension0021': land_cat_names,
                      'lat': (('west_east', 'south_north'),
                              dict_runs[WRFdomain]['CLAT'][0, ...].values),
                      'lon': (('west_east', 'south_north'),
                              dict_runs[WRFdomain]['CLONG'][0, ...].values)}
        dict_runs[WRFdomain] = dict_runs[WRFdomain].assign_coords(new_coords)
        dict_runs[WRFdomain] = dict_runs[WRFdomain].rename(
            {'z-dimension0021': 'PFT'})
    return(dict_runs)


if __name__ == '__main__':
    test_get_landuse_to_xarray = True
    test_get_data_file = False
    test_get_xarray = False
    test_get_xarray_daily = False
    test_merge = True

    if test_get_landuse_to_xarray:
        dict_runs = yatir_landuse_to_xarray()

    if test_get_data_file:
        roughness_data_dir = os.path.join('/', 'Users',
                                          'tim', 'work',
                                          'Data', 'SummenWRF',
                                          'yatir', 'roughness_experiments')
        fnames = {'ctl': 'energyfluxes_control_d03.nc',
                  'yatir': 'energyfluxes_yatir_run_d03.nc',
                  'vegparm_z0_10': 'energyfluxes_vegparm_z0_10_d03.nc',
                  'landuse_z0_500': 'energyfluxes_landuse_z0_500_d03.nc'}
        data_paths = {k: os.path.join(roughness_data_dir, v)
                      for k, v in fnames.items()}
        dsgv, dsxr = get_data_file(data_paths)

    if test_get_xarray:
        # LHctl = yatir_to_xarray(os.path.join('/', 'Users', 'tim', 'work',
        #                                      'Data', 'SummenWRF', 'yatir',
        #                                      'LH_d03_yatirZ50.nc'),
        #                         varname='LH',
        #                         groupname='ctl',
        #                         timerange=10)
        # LHytr = yatir_to_xarray(os.path.join('/', 'Users', 'tim', 'work',
        #                                      'Data', 'SummenWRF', 'yatir',
        #                                      'LH_d03_yatirZ50.nc'),
        #                         varname='LH',
        #                         groupname='yatirZ050',
        #                         timerange=10)

        ctl = yatir_WRF_to_xarray('/Users/tim/work/Data/SummenWRF/yatir/fluxes_ctl_run_d03.nc')
        ytr = yatir_WRF_to_xarray('/Users/tim/work/Data/SummenWRF/yatir/fluxes_yatir_run_d03.nc')

    if test_get_xarray_daily:
        ctl = yatir_WRF_to_xarray('/Users/tim/work/Data/SummenWRF/yatir/fluxes_ctl_run_d03.nc')
        ytr = yatir_WRF_to_xarray('/Users/tim/work/Data/SummenWRF/yatir/fluxes_yatir_run_d03.nc')
        ctlday = WRF_daily_daylight_avg('/Users/tim/work/Data/SummenWRF/yatir/fluxes_ctl_run_d03.nc')
        ytrday = WRF_daily_daylight_avg('/Users/tim/work/Data/SummenWRF/yatir/fluxes_yatir_run_d03.nc')
        dims = define_dims(ctl)

    if test_merge:
        ctlday, ytrday = merge_yatir_fluxes_landuse()
