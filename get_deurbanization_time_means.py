from datetime import datetime
import os
import numpy as np
import xarray as xr
import pandas as pd


deurbdir = "/global/cscratch1/sd/twhilton/WRFv4.0_Sensitivity/restored/deurb/"
ctldir = "/global/cscratch1/sd/twhilton/WRFv4.0_Sensitivity/restored/control/"
outdir = "/global/cscratch1/sd/twhilton/deurbanization_output_collected"


def get_Tmin_Tmax(ds):
    """calculate Tmin and Tmax from a dataset's T2 variable

    ARGS:
        ds: xarray dataset of WRF output; must contain dataarray T2
    """
    T2_daily = ds['T2'].groupby(group='doy')
    daily_Tmin = T2_daily.min(skipna=True)
    daily_Tmax = T2_daily.max(skipna=True)
    avg_Tmin = daily_Tmin.mean('doy')
    avg_Tmax = daily_Tmax.mean('doy')
    return(avg_Tmin, avg_Tmax)


def get_time_means(fname, run):
    ds = xr.open_dataset(fname, decode_times=False)
    for this_var in ds.data_vars:
        # WRF used varous values in the neighborgood of 9e36 for
        # missing values, but did not set the _FillValue attribute.
        # Replace these values with NaN now.
        try:
            print('replacing {}'.format(this_var))
            if this_var != "Times":
                ds[this_var] = ds[this_var].where(ds[this_var] < 1e33)
        except TypeError as e:
            print("error replacing missing values in " + this_var)
            raise e
    print('data read: ', datetime.now())
    timestamps = pd.to_datetime(ds.Times.astype(str).values.tolist(),
                                format='%Y-%m-%d_%H:%M:%S')
    ds = ds.assign_coords(month=xr.DataArray(timestamps.month.values,
                                             coords=[ds.Time],
                                             dims=['Time']),
                          doy=xr.DataArray(timestamps.dayofyear.values,
                                           coords=[ds.Time],
                                           dims=['Time']))
    print('starting groupby and mean: ', datetime.now())
    Tmin, Tmax = get_Tmin_Tmax(ds)
    month_grp = ds.groupby(group='month')
    month_means = month_grp.mean('Time', skipna=True)
    means = month_means.mean('month')
    means = means.assign_coords(longitude=ds['XLONG'][0, ...])
    means = means.assign_coords(latitude=ds['XLAT'][0, ...])
    means['Tmin'] = Tmin
    means['Tmax'] = Tmax
    means = means.drop(('doy', 'Time', 'month'))
    means.to_netcdf(path=os.path.join(outdir,
                                      'd02_{}_month_mean_vals_dbg.nc'.format(run)))
    print('done: ', datetime.now())
    return(means, ds)


if __name__ == "__main__":
    paths = {'ctl':os.path.join(outdir, 'ctl_d02.nc'),
             'deurb':os.path.join(outdir, 'deurb_d02.nc')}
    print('start: ', datetime.now())
    means = {run: get_time_means(fname, run)
             for run, fname in paths.items()}
