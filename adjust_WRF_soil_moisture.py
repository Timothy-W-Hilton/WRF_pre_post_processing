"""reduce soil moisture in met_em* netCDF4 files

met_em*.nc are the WRF driver data produced by met_grid.exe (of the
`WRF preprocessing system (WPS)<http://www2.mmm.ucar.edu/wrf/users/>`.
This module reduces the soil moisture fields in met_em* by an
arbitrary constant factor.

The motiviation behind this is testing WRF sensitivity to soil
moisture.

Timothy W. Hilton <thilton@ucmerced.edu>
12 Oct 2017
"""

import netCDF4
import numpy as np
import numpy.ma as ma
import glob
import os


def reduce_soil_moisture(fname, f, soil_moist_vars=None,
                         land_sea_var='LANDMASK'):
    """multiply WPS soil moisture values in a netcdf file by a constant factor

    WPS is the `WRF preprocessing system
    <http://www2.mmm.ucar.edu/wrf/users/>`.

    ARGS
    fname (character): full path to a met_em* file produced by the WPS.
    f (real): the constant factor to multiply into soil moisture values
    soil_moist_vars (list of strings): names of the soil moisture
       variables to multiply.  Default is ['SM100200', 'SM040100',
       'SM010040', 'SM000010']
    land_sea_var (character): name of the netCDF variable containing
       the land-sea mask for the dataset

    NOTES
    The netCDF file is altered in place, so the caller must have write
    priveleges for the file.
    """
    if soil_moist_vars is None:
        soil_moist_vars = ['SM100200', 'SM040100', 'SM010040', 'SM000010']

    nc = netCDF4.Dataset(fname, 'a')  # open in append mode
    land_sea = nc.variables[land_sea_var][...]
    # WRF land/sea mask is a float variable; water cells are 0.0, land
    # cells are 1.0.  Create a mask that is TRUE for water cells,
    # FALSE for land cells
    is_water = ma.masked_values(land_sea, 0.0).mask
    for this_var in soil_moist_vars:
        sm_adjusted = nc.variables[this_var][...] * f
        if sm_adjusted.ndim == 4:
            is_water = np.tile(is_water[:, np.newaxis, :, :], (1, 4, 1, 1))
        # set water cells back to 1.0
        sm_adjusted[is_water] = 1.0
        nc.variables[this_var][...] = sm_adjusted
        print "modified {this_var} in {fname}".format(this_var=this_var,
                                                      fname=fname)
    nc.history = ("created by metgrid_exe.  Soil moisture adjusted by"
                  " a factor of {f} by reduce_soil_moisture python"
                  " module.".format(f=f))
    nc.close()


def gather_files(dir, key):
    """return a list of all files in a directory that match a wildcard

    ARGS
    dir (character string): full path of directory to search
    key (character string): the search pattern

    RETURNS
    a list containing the full path of all files in dir whose names
    match key
    """
    return(glob.glob(os.path.join(dir, key)))


if __name__ == "__main__":
    wps_dir = os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton',
                           'WRFv3.9_Sensitivity',
                           'WRFv3.9_Sensitivity_DrySoil', 'WPS', '')
    wrf_run_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton', 'WRFv3.9_Sensitivity',
                               'met_em_tests', 'Dry',
                               'WRFV3', 'run')
    met_em_files = gather_files(wps_dir, "met_em*.nc")
    wrf_restart_files = gather_files(wrf_run_dir, "wrfrst_d*.nc")
    wrf_input_files = gather_files(wrf_run_dir, "wrfinput_d*")
    for this_file in wrf_input_files:
        reduce_soil_moisture(this_file, 0.05,
                             soil_moist_vars=("SMOIS", ))
