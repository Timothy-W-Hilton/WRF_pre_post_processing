import datetime
import os
import numpy.ma as ma
from plot_diff import var_diff, VarDiffPlotter
import means_diff_test

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
        vd.to_netcdf(os.path.join('/', 'global', 'cscratch1', 'sd',
                                  'twhilton',
                                  '{}_d{:02d}_RWurban.nc'.format(varname,
                                                                 DOMAIN)))
        print('done reading files ({})'.format(datetime.datetime.now() - t0))
        # vd.mask_land_or_water(mask_water=False)
    else:
        vd = var_diff(
            ncfile=os.path.join(
                cscratchdir,
                # '/Users/tim/work/Data/SummenWRF/',
                '{varname}_d{DOMAIN:02d}_RWurban.nc'.format(
                    varname=varname, DOMAIN=DOMAIN)))
        if vd.p is None:
            vd.get_significance_mask(significance=0.95, adj_autocorr=True)
    # vd = is_foggy_obrien_2013(vd)
    pfx = '7day_urban'
    # for this_series in ['all_tstamps', 'time_avg']:
    for this_series in ['time_avg']:
        if this_series == 'all_tstamps':
            t_end = 1
            pfx = 'redwoodsurban'
            time_title_str = None
        else:
            t_end = 1
            pfx = pfx + '_timeavg'
            vd.aggregate_time(time_avg=True)
            time_title_str = 'June 2009'
        for this_t in range(0, t_end):  #
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
                               mask=(vd.p < 0.95))
                               # vmin=-0.0000001,
                               # vmax=1.0000001)

    print('done driver ({})'.format(datetime.datetime.now() - t0))
