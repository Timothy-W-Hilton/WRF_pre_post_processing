import datetime
import os
import wrf
import netCDF4


def netCDF4_bigfile(wrf_dir):
    """timing test for reading a ~700 MB file using netCDF4.Dataset()
    """
    fname = 'HFX_d02_all.nc'
    t0 = datetime.datetime.now()
    print('starting {} read ({})'.format(fname, t0))
    bigfile = netCDF4.Dataset(os.path.join(wrf_dir, fname), 'r')
    # wrfv = wrf.getvar(bigfile, 'HFX', wrf.ALL_TIMES)
    netcdf4v = bigfile.variables['HFX'][...]
    bigfile.close()
    print('done {} ({})'.format(fname, datetime.datetime.now() - t0))
    return(netcdf4v)


def wrf_ncrcat_output(wrf_dir):
    """try to read the output of ncrcat using wrf.getvar()

    This works, but only if the ncrcat output includes the Times
    *variable*, in addition to the Time *dimension* (which is included
    by default)
    """
    fname = "HFX_1day_notime.nc"
    t0 = datetime.datetime.now()
    print('starting {} read with wrf.getvar() ({})'.format(fname, t0))
    bigfile = netCDF4.Dataset(os.path.join(wrf_dir, fname), 'r')
    wrfv = wrf.getvar(bigfile, 'HFX', wrf.ALL_TIMES)
    # netcdf4v = bigfile.variables['HFX'][...]
    bigfile.close()
    print('done {} ({})'.format(fname, datetime.datetime.now() - t0))
    return(wrfv)


if __name__ == "__main__":

    wrf_dir = os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton',
                           'WRFv4.0_Sensitivity', 'WRFCLMv4.0_NCEPDOEp2',
                           'WRFV4', 'run', 'summen_2005_ctl_NCEPDOE')
    wrfv = wrf_ncrcat_output(wrf_dir)
