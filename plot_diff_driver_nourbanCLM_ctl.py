import datetime
import os
import numpy as np
from plot_diff import var_diff, VarDiffPlotter
import sys

# ----- settings -----
DOMAIN = 2
PVAL = 0.05
READ_DATA = True
if len(sys.argv) > 1:
    varname = sys.argv[1]
else:
    # varname = 'HFX'
    varname = 'LH'
    # varname = 'LCL'
    # varname = 'fogpresent'
    # varname = 'SMOIS'
    # varname = 'fogbase'
    # varname = 'fogpct'
    # varname = 'QCLOUD'
    # varname = 'uvmet10'
    # varname = 'wa'
    # varname = 'ctt'  # cloud top temperature
    # varname = 'LU_INDEX'
print('varname: {}'.format(varname))

if __name__ == "__main__":
    t0 = datetime.datetime.now()

    cscratchdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton')
    rootdir = os.path.join(cscratchdir, 'WRFv4.0_Sensitivity')
    ctl_dir = os.path.join(rootdir, 'WRFCLMv4.0_NCEPDOEp2', 'WRFV4',
                           'run', 'summen_2005_ctl_NCEPDOE')
    nourb_dir = os.path.join(rootdir,
                             'WRFCLMv4.0_NCEPDOEp2_deurbanized', 'WRFV4',
                             'run', 'summen_2005_deurbanized_NCEPDOE')
    out_dir = os.path.join(cscratchdir, 'plots_temporary')

    if READ_DATA:
        # wildcard_pat = "*d{:02d}02_2009-06-\{0[0-9],1[01]\}".format(DOMAIN)
        wildcard_pat = "*d{:02d}_2005-*".format(DOMAIN)
        re_pat = "d{:02d}_all\\.nc".format(DOMAIN)
        vd = var_diff(os.path.join(ctl_dir, re_pat),
                      os.path.join(nourb_dir, re_pat),
                      label_A='ctl',
                      label_B='no_urban_CLM',
                      varname=varname)
        vd.read_files()
        # for k in vd.data.keys():
        #     z_ax = vd.var_axes.index('Lay')
        #     idx = [slice(None), ] * vd.data[k].ndim
        #     idx[z_ax] = slice(0, 10)
        #     vd.data[k] = np.mean(vd.data[k][idx], axis=z_ax)
        # vd.var_axes.pop(vd.var_axes.index('Lay'))
        vd.get_significance_mask(significance=0.95, adj_autocorr=True)
        vd.to_netcdf(os.path.join('/', 'global', 'cscratch1', 'sd',
                                  'twhilton',
                                  '{}_d{:02d}_CLM_nourban.nc'.format(
                                      varname,
                                      DOMAIN)))
        print('done reading files ({})'.format(datetime.datetime.now() - t0))
        # vd.mask_land_or_water(mask_water=False)
    else:
        vd = var_diff(
            ncfile=os.path.join(
                cscratchdir,
                # '/Users/tim/work/Data/SummenWRF/',
                '{varname}_d{DOMAIN:02d}_CLM_nourban.nc'.format(
                    varname=varname, DOMAIN=DOMAIN)))
        # for k in vd.data.keys():
        #     vd.data[k] = vd.data[k] * 100.0
        if vd.p is None:
            vd.get_significance_mask(significance=0.95, adj_autocorr=True)
    # vd = is_foggy_obrien_2013(vd)
    pfx = 'nourbanCLM_{:0.0f}CI'.format((1.0 - PVAL) * 100.0)
    # for this_series in ['all_tstamps', 'time_avg']:
    for this_series in ['time_avg']:
        if this_series == 'all_tstamps':
            t_end = 1
            pfx = 'nourbanCLM_1day'
            time_title_str = None
        else:
            t_end = 1
            pfx = pfx + '_timeavg'
            vd.aggregate_time(time_avg=True)
            time_title_str = ''
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
            fig = plotter.plot(cb_orientation=cb_orientation,
                               vmin=None,
                               vmax=None,
                               mask=vd.p > PVAL)
                               # mask=None)


    print('done driver ({})'.format(datetime.datetime.now() - t0))
