import datetime
import os
import xarray as xr
import re
import wrf
import netCDF4

n_files = 200

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
    ds = xr.open_mfdataset(matching_files[0])
    print('done opening first file'.format(
        datetime.datetime.now() - t0))

    t0 = datetime.datetime.now()
    all_vars = list(ds.data_vars.keys())
    ds.close()
    keep_vars = ['Times', 'HFX', 'LH', 'QCLOUD']
    drop_vars = set(all_vars) - set(keep_vars)
    print('starting xarray.open_mfdataset read (with drop) ({})'.format(t0))
    ds = xr.open_mfdataset(matching_files[:n_files],
                           drop_variables=drop_vars)
    print('done xarray.open_mfdataset ({}) ({})'.format(
        n_files,
        datetime.datetime.now() - t0))

    # # this doesn't work
    # fname = 'HFX_d02_all.nc'
    # fname = 'foo1day.nc'
    # fname = 'metem_2005_ctl_NCEPDOE_d02_2005-06-01_00:00:00.nc'
    # print('starting {} read ({})'.format(fname, t0))
    # bigfile = netCDF4.Dataset(os.path.join(wrf_dir, fname), 'r')
    # wrfv = wrf.getvar(bigfile, 'HFX', wrf.ALL_TIMES)
    # # netcdf4v = bigfile.variables['HFX'][...]
    # bigfile.close()
    # print('done {} ({})'.format(fname, datetime.datetime.now() - t0))
