import datetime
import os
import numpy.ma as ma
from plot_diff import var_diff, VarDiffPlotter
import means_diff_test

from map_tools_twh.map_tools_twh import Yatir_WRF_prj
from map_tools_twh.map_tools_twh import Yatir_WRF_Mapper

DOMAIN = 2

if __name__ == "__main__":
    t0 = datetime.datetime.now()

    cscratchdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton')
    rootdir = os.path.join(cscratchdir, 'WRFv4.1_Experiments')
    ctl_dir = os.path.join(rootdir,
                           ('WRFv4.1_yatir_2015Aug_CTL_'
                            'NCEPFNL_NOAHMP_ndom3'),
                           'WRFV4', 'run',
                           'yatir_2015Aug_CTL_NCEPFNL_NOAHMP_ndom3')
    exp_dir = os.path.join(rootdir,
                           ('WRFv4.1_yatir_2015Aug_yatirparams_'
                            'Z050_NCEPFNL_NOAHMP_ndom3'),
                           'WRFV4', 'run',
                           'yatir_2015Aug_yatirparams_NCEPFNL_NOAHMP_ndom3')
    out_dir = os.path.join(cscratchdir, 'plots_temporary')

    # insignificant_pixels = ma.masked_invalid(means_diff_test.main(0.95)).mask

    varname = 'HFX'

    read_data = False
    if read_data:
        wildcard_pat = "*d{:02d}_2015-08-1[567]*".format(DOMAIN)
        re_pat = "d{:02d}_2015-08-1[567].*\\.nc".format(DOMAIN)
        vd = var_diff(os.path.join(ctl_dir, re_pat),
                      os.path.join(exp_dir, re_pat),
                      label_A='ctl',
                      label_B='yatirZ050',
                      varname=varname)
        vd.read_files()
        vd.get_significance_mask(significance=0.95, adj_autocorr=True)
        vd.to_netcdf(os.path.join('/', 'global', 'cscratch1', 'sd',
                                  'twhilton',
                                  '{}_d{:02d}_yatirZ50.nc'.format(varname,
                                                                 DOMAIN)))
        print('done reading files ({})'.format(datetime.datetime.now() - t0))
        # vd.mask_land_or_water(mask_water=False)
    else:
        vd = var_diff(
            ncfile=os.path.join(
                cscratchdir,
                # '/Users/tim/work/Data/SummenWRF/',
                '{varname}_d{DOMAIN:02d}_yatirZ50.nc'.format(
                    varname=varname, DOMAIN=DOMAIN)))
        if vd.p is None:
            vd.get_significance_mask(significance=0.95, adj_autocorr=True)
    # vd = is_foggy_obrien_2013(vd)
    pfx = 'yatirZ050'
    # for this_series in ['all_tstamps', 'time_avg']:
    for this_series in ['time_avg']:
        if this_series == 'all_tstamps':
            t_end = 1
            pfx = 'yatirZ050'
            time_title_str = None
        else:
            t_end = 1
            pfx = pfx + '_timeavg'
            vd.aggregate_time(time_avg=True)
            time_title_str = '15-17 Aug 2015'
        for this_t in range(0, t_end):  #
            plotter = VarDiffPlotter(vd, t_idx=this_t, layer=0,
                                     domain=DOMAIN,
                                     pfx=pfx,
                                     savedir=out_dir,
                                     time_title_str=time_title_str,
                                     map_prj=Yatir_WRF_prj(),
                                     mapper=Yatir_WRF_Mapper)
            if DOMAIN == 1:
                cb_orientation = 'horizontal'
            else:
                cb_orientation = 'vertical'
            fig = plotter.plot(cb_orientation=cb_orientation,
                               mask=(vd.p < 0.95))
                               # vmin=-0.0000001,
                               # vmax=1.0000001)

    print('done driver ({})'.format(datetime.datetime.now() - t0))
