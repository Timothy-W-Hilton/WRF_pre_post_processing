from datetime import datetime
import os
import numpy as np
import xarray as xr
import pandas as pd

def get_time_means(fname, run):
    ds = xr.open_dataset(fname)
    print('data read: ', datetime.now())
    timestamps = pd.to_datetime(ds.Times.astype(str).values.tolist(),
                                format='%Y-%m-%d_%H:%M:%S')
    ds = ds.assign_coords(month=xr.DataArray(timestamps.month.values,
                                             coords=[ds.Time],
                                             dims=['Time']))
    print('starting groupby and mean: ', datetime.now())
    month_grp = ds.groupby(group='month')
    month_means = month_grp.mean('Time')
    means = month_means.mean('month')
    means = means.assign_coords(longitude=ds['XLONG'][0, ...])
    means = means.assign_coords(latitude=ds['XLAT'][0, ...])
    means.to_netcdf(path=os.path.join(os.path.dirname(fname),
                                      'd02_{}_10day_mean_vals.nc'.format(run)))
    print('done: ', datetime.now())
    return(means)

if __name__ == "__main__":
    paths = {'ctl':os.path.join('/', 'global', 'cscratch1', 'sd',
                                'twhilton', 'WRFv4.0_Sensitivity',
                                'WRFCLMv4.0_NCEPDOEp2', 'WRFV4', 'run',
                                'summen_2005_ctl_NCEPDOE', 'd02_all.nc'),
             'deurb':os.path.join('/', 'global', 'cscratch1', 'sd',
                                  'twhilton', 'WRFv4.0_Sensitivity',
                                  'WRFCLMv4.0_NCEPDOEp2_deurbanized',
                                  'WRFV4', 'run',
                                  'summen_2005_deurbanized_NCEPDOE',
                                  'd02_all.nc')}
    paths_10day = {'ctl': os.path.join('/', 'global', 'cscratch1',
                                       'sd', 'twhilton',
                                       'deurbanization_output_collected',
                                       'ctl_d02_10day.nc'),
                   'deurb':os.path.join('/', 'global', 'cscratch1',
                                        'sd', 'twhilton',
                                        'deurbanization_output_collected',
                                        'deurb_d02_10day.nc')}
    print('start: ', datetime.now())
    means = {run: get_time_means(fname, run)
             for run, fname in paths_10day.items()}
