import os
import xarray as xr
import netCDF4
import wrf


def create_merged_postprocessed_netcdf(collected_data_directory,
                                       allvars_fname,
                                       vars_derived,
                                       vars_raw,
                                       run_key):

    fname = os.path.join(collected_data_directory,
                         allvars_fname[run_key])
    nc = netCDF4.Dataset(filename=fname, mode='r')

    for k in vars_derived:
        try:
            print('calculating {}'.format(k))
            wrf.getvar(wrfin=nc,
                       varname=k,
                       timeidx=wrf.ALL_TIMES)
        except ValueError as e:
            print('error calculating {}'.format(k))
            raise(e)

    ds_derived = xr.Dataset({k: wrf.getvar(wrfin=nc,
                                           varname=k,
                                           timeidx=wrf.ALL_TIMES)
                             for k in vars_derived})

    nc.close()

    # format variable metadata for compatibility with
    # xarray.Dataset.to_netcdf().
    for this_var in ds_derived.data_vars:
        for k, v in ds_derived[this_var].attrs.items():
            if type(v) is wrf.projection.LambertConformal:
                # cast projection information into string -
                # wrf.projection objects throw an error in to_netcdf()
                ds_derived[this_var].attrs[k] = str(v)
            elif k == 'coordinates':
                # (2) rename variables' "coordinates" attribute to
                # 'coordinates_attr'. See
                # https://github.com/pydata/xarray/issues/1809)
                ds_derived[this_var].attrs[
                    'coordinates_attr'] = ds_derived[this_var].attrs.pop(k)

    ds_raw = xr.open_dataset(fname)
    # set(vars_raw)) list ds_raw first and set compat to override
    # because getvar() replaces string timestamps with numerical
    # timestamps, and merge() then sees that ds_raw has strings in
    # XTIME and ds_derived has numbers and decides they are
    # incompatible.  setting compat to override keeps the values from
    # the first dataset (the strings in raw).
    ds_combined = xr.merge([ds_raw[vars_raw], ds_derived], compat='override')
    ds_combined.to_netcdf(path=os.path.join(
        os.path.dirname(fname),
        '{}_d03_VWCx2_postprocessed.nc'.format(run_key)))


if __name__ == "__main__":

    collected_data_directory = os.path.join('/', 'global',
                                            'cscratch1', 'sd', 'twhilton',
                                            'yatir_output_collected',
                                            'wetsoil')
    allvars_fname = {'ytr': 'ytr_d03_allvars.nc',
                     'ctl': 'ctl_d03_allvars.nc'}
    allvars_fname = {'ytr': 'yatir_run_d03_diag_TP_VWCx2.nc',
                     'ctl': 'yatir_run_d03_diag_TP_VWCx1.nc'}

    vars_derived = ['height', 'height_agl', 'ter', 'zstag', 'uvmet_wspd_wdir',
                    'uvmet10_wspd_wdir', 'theta']
    vars_raw = ['Times', 'XLONG', 'XLAT', 'HFX', 'LH', 'FIRA', 'FSA',
                'GPP', 'SWDOWN', 'FVEG', 'RSSUN', 'RSSHA', 'W', 'T', 'P']

    for k in ['ytr', 'ctl']:
        create_merged_postprocessed_netcdf(collected_data_directory,
                                           allvars_fname,
                                           vars_derived,
                                           vars_raw,
                                           k)
