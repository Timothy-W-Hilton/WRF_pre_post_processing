import datetime
import numpy as np
import os
import xarray as xr
import copy
from plot_diff import var_diff, VarDiffPlotter

DOMAIN = 2

if __name__ == "__main__":
    t0 = datetime.datetime.now()

    cscratchdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton')
    rootdir = os.path.join(cscratchdir, 'WRFv3.9_Sensitivity')
    ctl_dir = os.path.join(rootdir, 'WRFV3_Ctl',
                           'run', 'summen_sensitivity_ctl')
    urb_dir = os.path.join(rootdir, 'WRFV3_redwoodsurban',
                           'run', 'summen_redwoodsurban')
    out_dir = os.path.join(cscratchdir, 'plots_temporary')

    # insignificant_pixels = ma.masked_invalid(means_diff_test.main(0.95)).mask

    varname = 'fogpresent'

    read_data = False
    if read_data:
        # wildcard_pat = "*d{:02d}02_2009-06-\{0[0-9],1[01]\}".format(DOMAIN)
        wildcard_pat = "*d{:02d}_2009-06-*".format(DOMAIN)
        vd = var_diff(os.path.join(ctl_dir, wildcard_pat),
                      os.path.join(urb_dir, wildcard_pat),
                      label_A='ctl',
                      label_B='urbanRW',
                      varname=varname)
        vd.read_files()
        vd.get_significance_mask(significance=0.95, adj_autocorr=True)
        # vd.to_netcdf(os.path.join('/', 'global', 'cscratch1', 'sd',
        #                           'twhilton',
        #                           '{}_d{:02d}_RWurban.nc'.format(varname,
        #                                                          DOMAIN)))
        print('done reading files ({})'.format(datetime.datetime.now() - t0))
        # vd.mask_land_or_water(mask_water=False)
    else:
        vd = var_diff(
            ncfile=os.path.join(
                # cscratchdir,
                '/Users/tim/work/Data/SummenWRF/',
                '{varname}_d{DOMAIN:02d}_RWurban.nc'.format(
                    varname=varname, DOMAIN=DOMAIN)))
        if vd.p is None:
            vd.get_significance_mask(significance=0.95, adj_autocorr=True)
    # vd = is_foggy_obrien_2013(vd)

    if False:
        vd.get_pval_timeseries(interval_hrs=336.0)
        # now need to write pval timeseries to netCDF or something....
        pvals_xr = xr.Dataset({'p': (('t', 'x', 'y'),
                                     vd.pvals_series),
                               'mean_diff': (('t', 'x', 'y'),
                                             vd.mean_diff_series),
                               'idx': (('t', ), vd.idx_pvals_series),
                               'Latitude': (('x', 'y'),
                                            vd.lat),
                               'Longitude': (('x', 'y'),
                                             vd.lon)},
                              coords={'t': vd.t_pvals_series,
                                      'x': range(vd.pvals_series.shape[1]),
                                      'y': range(vd.pvals_series.shape[2])})
        ncfname = '/Users/tim/work/Data/SummenWRF/pvals_series.nc'
        os.remove(ncfname)
        pvals_xr.to_netcdf(ncfname)
    else:
        pvals_xr = xr.open_dataset(os.path.join('/', 'Users', 'tim',
                                                'work', 'Data', 'SummenWRF',
                                                'pvals_series.nc'))

    data_orig = copy.copy(vd.data)
    for i, this_t in enumerate(pvals_xr.t.values):
        vd.d = pvals_xr.mean_diff.data[i, ...]
        vd.abs_max = np.nanmax(np.abs(pvals_xr.mean_diff.data))

        for k in data_orig.keys():
            vd.data[k] = np.mean(data_orig[k][:pvals_xr.idx.values[i], ...],
                                 axis=0,
                                 keepdims=True)
        plotter = VarDiffPlotter(vd,
                                 fig_type='png',
                                 t_idx=0,
                                 layer=0,
                                 domain=DOMAIN,
                                 pfx='pval_test',
                                 savedir='/Users/tim/work/Plots/Summen/',
                                 time_title_str=np.datetime_as_string(
                                     this_t, unit='m'),
                                 show_title=False)
        fig = plotter.plot(cb_orientation='vertical',
                           mask=pvals_xr.p.data[i, ...] < 0.99)
