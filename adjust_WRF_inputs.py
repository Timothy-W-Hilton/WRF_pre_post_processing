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


def make_all_land_urban(fname):
    """set all land use in a real.exe wrfinput netcdf file to urban

    real.exe takes meteorological data prepared by WPS and creates
    wrfinput files.  WPS is the `WRF preprocessing system
    <http://www2.mmm.ucar.edu/wrf/users/>`.

    ARGS
    fname (character): full path to a met_em* file produced by the WPS.

    NOTES
    - The netCDF file is altered in place, so the caller must have
      write priveleges for the file.
    - There are several land use datasets available to WRF.  The
      default is MODIS.  See Ch. 3 of the WRF user guide, section
      "Selecting Between USGS and MODIS-based Land Use
      Classifications".  See also file LANDUSE.TBL in the WRF run
      directory.
    - Currently doesn't check the landmask data in the netCDF file.
      It assumes that WPS correctly sets LU_INDEX for all water pixels
      to the value contained in the ISWATER global netCDF attribute.
      The value of ISWATER depends on the land use dataset selected.
    """
    nc = netCDF4.Dataset(fname, 'a')  # open in append mode
    # set LU_INDEX to urban for all non-water pixels.  This approach
    #
    landuse = nc.variables['LU_INDEX'][...]
    landuse[np.where(landuse != nc.ISWATER)] = nc.ISURBAN
    nc.variables['LU_INDEX'][...] = landuse
    print "modified LU_INDEX in {fname}".format(fname=fname)
    # LANDUSEF specifies a fraction for all land use types in each
    # pixel.  Set all urban fraction in LANDUSEF to (1.0 -
    # LANDUSEF[LU_water]) and all other non-water fractions to 0.0.
    all_land_use = np.arange(1, nc.NUM_LAND_CAT)
    non_water_non_urban = np.setdiff1d(all_land_use, [nc.ISWATER, nc.ISURBAN])
    # need to subtract 1 from land use codes to translante them to
    # zero-based array indices
    nc.variables['LANDUSEF'][:, (non_water_non_urban - 1), ...] = 0.0
    non_water_fraction = (
        1.0 - nc.variables['LANDUSEF'][:, (nc.ISWATER - 1), ...])
    nc.variables['LANDUSEF'][:, (nc.ISURBAN - 1), ...] = non_water_fraction
    # TODO change module name to reflect that it
    # now can alter landuse
    nc.history = ("created by metgrid_exe.  All non-waterland use set "
                  "to urban by adjust_WRF_soil_moisture python module.")
    nc.close()


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
