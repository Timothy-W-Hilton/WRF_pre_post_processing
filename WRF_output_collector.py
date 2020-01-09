import os
import pandas as pd
from nco import Nco
import timeit


def get_fnames_list(prefix, t0, t_end):

    # T is the alias for minutes(!)
    tstamps = pd.date_range(start=t0, end=t_end, freq='30T')
    return(list(tstamps.strftime(prefix +
                                 '%Y-%m-%d_%H:%M:%S.nc')))


def concatenate_files(dirpath, files_list, fname_out, var_list=None):
    """use ncrcat to concatenate WRF files into one large file

    This is useful for using, for example, wrf.getvar() to calculate
    diagnostics.

    """
    if var_list is None:
        opts = None
    else:
        opts = ['-v ' + var_list]
    start = timeit.default_timer()
    files_full_paths = [os.path.join(dirpath, f) for f in files_list]
    print('starting ' + os.path.basename(fname_out))
    Nco().ncrcat(input=files_full_paths,
                 output=fname_out,
                 options=opts)
    end = timeit.default_timer()
    print('concatenated dataset ({:.0f} s)'.format(end-start))


if __name__ == "__main__":
    dirs = {'ytr': os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton',
                                'WRFv4.1_Experiments',
                                ('WRFv4.1_yatir_2015Aug_yatirparams_'
                                 'Z050_NCEPFNL_NOAHMP_ndom3_diagnostic'),
                                'WRFV4', 'run',
                                ('yatir_2015Aug_yatirparams_'
                                 'NCEPFNL_NOAHMP_ndom3')),
            'ctl': os.path.join('/', 'global', 'cscratch1', 'sd',
                                'twhilton', 'WRFv4.1_Experiments',
                                ('WRFv4.1_yatir_2015Aug_CTL_'
                                 'NCEPFNL_NOAHMP_ndom3'),
                                'WRFV4', 'run',
                                'yatir_2015Aug_CTL_NCEPFNL_NOAHMP_ndom3')}
    fname_prefix = {'ytr': ('metem_yatir_2015Aug_yatirparams_NCEPFNL_'
                            'NOAHMP_ndom3_d03_'),
                    'ctl': ('metem_yatir_2015Aug_CTL_'
                            'NCEPFNL_NOAHMP_ndom3_d03_')}

    t0 = pd.datetime(2015, 8, 15, 0, 0, 0)
    t_end = pd.datetime(2015, 8, 15, 2, 0, 0)
    fnames_list = {k: get_fnames_list(v, t0, t_end) for k, v in
                   fname_prefix.items()}
    vars = "Times,XLONG,XLAT,W,T,P,T2V,T2B,TAH,TV,TH2,T2"

    for k in dirs.keys():
        concatenate_files(dirs[k],
                          fnames_list[k],
                          os.path.join('/', 'global', 'cscratch1', 'sd',
                                       'twhilton', 'yatir_output_collected',
                                       '{}_d03_vartest.nc'.format(k)),
                          var_list=vars)
