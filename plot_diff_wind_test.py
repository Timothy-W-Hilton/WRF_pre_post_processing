import datetime
import os
import numpy as np
from plot_diff import var_diff, VarDiffPlotter

DOMAIN = 2

if __name__ == "__main__":
    t0 = datetime.datetime.now()

    cscratchdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton')
    rootdir = os.path.join(cscratchdir, 'WRFv4.0_Sensitivity')
    ctl_dir = os.path.join(rootdir, 'WRFCLMv4.0_NCEPDOEp2', 'WRFV4',
                           'run', 'summen_2005_ctl_NCEPDOE')
    dry_dir = os.path.join(rootdir,
                           'WRFCLMv4.0_NCEPDOEp2_drysoil', 'WRFV4',
                           'run', 'summen_2005_NCEPDOE_drysoil')

    out_dir = os.path.join(cscratchdir, 'plots_temporary')

    varname = 'uvmet'
    varname = 'uvmet10'
    # varname = 'wa'
    # varname = 'ctt'  # cloud top temperature

    read_data = False
    if read_data:
        # wildcard_pat = "*d{:02d}02_2009-06-\{0[0-9],1[01]\}".format(DOMAIN)
        wildcard_pat = "*d{:02d}_2005-06-0*".format(DOMAIN)
        vd = var_diff(os.path.join(ctl_dir, wildcard_pat),
                      os.path.join(dry_dir, wildcard_pat),
                      label_A='ctl',
                      label_B='drysoil_CLM',
                      varname=varname)
        vd.read_files()
        # for k in vd.data.keys():
        #     z_ax = vd.var_axes.index('Lay')
        #     idx = [slice(None), ] * vd.data[k].ndim
        #     idx[z_ax] = slice(0, 10)
        #     vd.data[k] = np.mean(vd.data[k][idx], axis=z_ax)
        # vd.var_axes.pop(vd.var_axes.index('Lay'))
        # vd.get_significance_mask(significance=0.95, adj_autocorr=True)
        vd.to_netcdf(os.path.join('/', 'global', 'cscratch1', 'sd',
                                  'twhilton',
                                  '{}_d{:02d}_CLM_drysoil.nc'.format(varname,
                                                                     DOMAIN)))
        print('done reading files ({})'.format(datetime.datetime.now() - t0))
    else:
        vd = var_diff(
            ncfile=os.path.join(
                '/', 'global', 'cscratch1', 'sd',
                'twhilton',
                '{varname}_d{DOMAIN:02d}_CLM_drysoil.nc'.format(
                    varname=varname, DOMAIN=DOMAIN)))
    # vd = is_foggy_obrien_2013(vd)
    # vd.get_significance_mask(significance=0.95, adj_autocorr=True)

    pfx = 'TEST_drysoilCLM'
    # for this_series in ['all_tstamps', 'time_avg']:
    for this_series in ['single']:
        if this_series == 'all_tstamps':
            t_end = 1
            pfx = 'drysoilCLM'
            time_title_str = None
        elif this_series == 'single':
            t_end = 200
        else:
            t_end = 1
            pfx = pfx + '_timeavg'
            vd.aggregate_time(time_avg=True)
            if varname is "fogpresent":
                print("multiplying fogpresent avg by 100 to get a percent")
                for k in vd.data.keys():
                    vd.data[k] = vd.data[k] * 100
            time_title_str = 'June 2005'
        for this_t in range(10, t_end):  #
            if this_series == 'single':
                time_title_str = vd.time[this_t].strftime('%Y-%m-%d_%H%M')
            plotter = VarDiffPlotter(vd, t_idx=this_t, layer=0,
                                     domain=DOMAIN,
                                     pfx=pfx,
                                     savedir=out_dir,
                                     time_title_str=time_title_str)
            if DOMAIN == 1:
                cb_orientation = 'horizontal'
            else:
                cb_orientation = 'vertical'
            fig = plotter.plot(cb_orientation=cb_orientation,
                               vmin=None,
                               vmax=None,
                               mask=None)

    print('done driver ({})'.format(datetime.datetime.now() - t0))
