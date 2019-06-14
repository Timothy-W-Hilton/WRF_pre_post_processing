import datetime
import os
import xarray as xr
import re
import wrf
import netCDF4

if __name__ == "__main__":

    wrf_dir = os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton',
                           'WRFv4.0_Sensitivity', 'WRFCLMv4.0_NCEPDOEp2',
                           'WRFV4', 'run', 'summen_2005_ctl_NCEPDOE')
    re_pat = "d02_2005-"

    matching_files = [f for f in filter(re.compile(re_pat).search,
                                        os.listdir(wrf_dir))]
    # put the full directory path back in
    matching_files = [os.path.join(wrf_dir, f)
                      for f in matching_files]
    matching_files = sorted(matching_files)

    t0 = datetime.datetime.now()
    print('starting xarray.open_mfdataset read ({})'.format(t0))
    ds = xr.open_mfdataset(matching_files)
    print('done xarray.open_mfdataset ({})'.format(
        datetime.datetime.now() - t0))

    t0 = datetime.datetime.now()
    print('starting HFX read ({})'.format(t0))
    hfx = ds['HFX'].values
    print('done HFX read ({})'.format(
        datetime.datetime.now() - t0))

    # this doesn't work
    fname = 'HFX_d02_all.nc'
    fname = 'foo1day.nc'
    fname = 'metem_2005_ctl_NCEPDOE_d02_2005-06-01_00:00:00.nc'
    print('starting {} read ({})'.format(fname, t0))
    bigfile = netCDF4.Dataset(os.path.join(wrf_dir, fname), 'r')
    wrfv = wrf.getvar(bigfile, 'HFX', wrf.ALL_TIMES)
    # netcdf4v = bigfile.variables['HFX'][...]
    bigfile.close()
    print('done {} ({})'.format(fname, datetime.datetime.now() - t0))
