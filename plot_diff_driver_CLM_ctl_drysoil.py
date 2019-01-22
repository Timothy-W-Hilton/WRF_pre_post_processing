import datetime
import os
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

    # varname='HFX'
    # varname='LH'
    # varname='LCL'
    # varname='fogpresent'
    # varname='SMOIS'
    # varname='fogbase'
    varname='fogpct'
    # varname='QCLOUD'
    # varname='uvmet'
    # varname='wa'
    # varname='LU_INDEX'

    read_data = True
    if read_data:
        # wildcard_pat = "*d{:02d}02_2009-06-\{0[0-9],1[01]\}".format(DOMAIN)
        wildcard_pat = "*d{:02d}_2005-*".format(DOMAIN)
        vd = var_diff(os.path.join(ctl_dir, wildcard_pat),
                      os.path.join(dry_dir, wildcard_pat),
                      label_A='ctl',
                      label_B='drysoil_CLM',
                      varname=varname)
        vd.read_files()
        vd.to_netcdf(os.path.join('/', 'global', 'cscratch1', 'sd',
                                  'twhilton',
                                  '{}_d{:02d}_CLM_drysoil.nc'.format(varname,
                                                                     DOMAIN)))
        print('done reading files ({})'.format(datetime.datetime.now() - t0))
        # vd.mask_land_or_water(mask_water=False)
    else:
        vd = var_diff(
            ncfile=os.path.join('/', 'global', 'cscratch1', 'sd',
                                'twhilton',
                                '{varname}_d{:02d}_CLM_drysoil.nc'.format(
                                    varname, DOMAIN)))
    # vd = is_foggy_obrien_2013(vd)
    pfx = 'drysoilCLM'
    # for this_series in ['all_tstamps', 'time_avg']:
    for this_series in ['time_avg']:
        if this_series == 'all_tstamps':
            t_end = 1
            pfx = 'drysoilCLM'
            time_title_str = None
        else:
            t_end = 1
            pfx = pfx + '_timeavg'
            vd.aggregate_time(time_avg=True)
            time_title_str = 'June 2005'
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
                               vmin=1,
                               vmax=21)

    print('done driver ({})'.format(datetime.datetime.now() - t0))
