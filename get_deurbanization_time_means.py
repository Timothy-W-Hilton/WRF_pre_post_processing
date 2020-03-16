from datetime import datetime
import os
import numpy as np
import xarray as xr
import pandas as pd


deurbdir = "/global/cscratch1/sd/twhilton/WRFv4.0_Sensitivity/restored/deurb/"
ctldir = "/global/cscratch1/sd/twhilton/WRFv4.0_Sensitivity/restored/control/"
outdir = "/global/cscratch1/sd/twhilton/deurbanization_output_collected"



def get_time_means(fname, run):
    ds = xr.open_dataset(fname, decode_times=False)
    print('data read: ', datetime.now())
    timestamps = pd.to_datetime(ds.Times.astype(str).values.tolist(),
                                format='%Y-%m-%d_%H:%M:%S')
    ds = ds.assign_coords(month=xr.DataArray(timestamps.month.values,
                                             coords=[ds.Time],
                                             dims=['Time']))
    print('starting groupby and mean: ', datetime.now())
    month_grp = ds.groupby(group='month')
    month_means = month_grp.mean('Time', skipna=True)
    means = month_means.mean('month')
    means = means.assign_coords(longitude=ds['XLONG'][0, ...])
    means = means.assign_coords(latitude=ds['XLAT'][0, ...])
    means.to_netcdf(path=os.path.join(outdir,
                                      'd02_{}_month_mean_vals.nc'.format(run)))
    print('done: ', datetime.now())
    return(means)


if __name__ == "__main__":
    paths = {'ctl':os.path.join(outdir, 'ctl_d02.nc'),
             'deurb':os.path.join(outdir, 'deurb_d02.nc')}
    # paths_10day = {'ctl': os.path.join('/', 'global', 'cscratch1',
    #                                    'sd', 'twhilton',
    #                                    'deurbanization_output_collected',
    #                                    'ctl_d02_10day.nc'),
    #                'deurb':os.path.join('/', 'global', 'cscratch1',
    #                                     'sd', 'twhilton',
    #                                     'deurbanization_output_collected',
    #                                     'deurb_d02_10day.nc')}
    print('start: ', datetime.now())
    means = {run: get_time_means(fname, run)
             for run, fname in paths.items()}
