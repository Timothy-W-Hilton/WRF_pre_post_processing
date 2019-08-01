"""reduce soil moisture in met_em* netCDF4 files

met_em*.nc are the WRF driver data produced by met_grid.exe (of the
`WRF preprocessing system (WPS)<http://www2.mmm.ucar.edu/wrf/users/>`.
This module reduces the soil moisture fields in met_em* by an
arbitrary constant factor.

The motiviation behind this is testing WRF sensitivity to soil
moisture.

Timothy W. Hilton <twhilton@ucsc.edu>
12 Oct 2017
"""

import netCDF4
import numpy as np
import numpy.ma as ma
import glob
import os
import geopy.distance

def deal_with_100pct_urban(landusef, f_urban):
    """replace urban land use in cells that are 100% urban

    These require a decision, because these cells have no non-urban
    land use to scale up.  The decision I am making is to find the
    nearest neighbor to the east with some non-urban land use, and use
    that cell.  I am using the nearest neighbor to the East because
    all of the 100% urban cells in California are in LA near the
    coast, so moving inland to search for a replacement seems
    reasonable.  I'm aware moving inland could introduce topography
    effects; hopefully these are minor.

    Another possible approach is to use the outer domain cell that
    contains the current inner domain cell.  That's a more complex
    calculation, so punting for now.
    """
    idx_all_urban = zip(*np.where(np.isclose(f_urban, 1.0)))
    for t, y, x in idx_all_urban:
        found_replacement = False
        xeast = x
        while not(found_replacement):
            xeast = xeast + 1  # try the next cell to the east by
                               # incrementing x index
            if (f_urban[0, y, xeast] < 0.99):
                found_replacement = True
        landusef[0, :, y, x] = landusef[0, :, y, xeast]
        f_urban[0, y, x] = f_urban[0, y, xeast]
        print("    replaced {} land use with {}".format((y, x), (y, xeast)))
    return(landusef, f_urban)


def remove_urban(fname_wrf):
    """remove urban landuse from WRF input files

    Set urban landuse fraction to 0.0, increase non-water land uses
    such that the sum of all land uses remains 1.0 and the relative
    proportions of non-water landuses are unchanged, and reset the
    land-use index to the highest-fraction remaining landuse in teh
    gridcell.
    """
    nc = netCDF4.Dataset(fname_wrf, 'a')  # open in append mode
    landusef = np.nan_to_num(np.copy(nc.variables['LANDUSEF'][...]))
    # this algorithm will fail for grid cells where f_urban is 1.0,
    # because the denominator in the adjustment factor will be zero.
    # Conceptually, these cells will need user intervention to decide
    # what to do, because there are no non-urban land covers to scale
    # up.
    f_urban = np.copy(landusef[:, nc.ISURBAN - 1, ...])
    if (np.isclose(f_urban, 1.0).any()):
        print(('{ncells} cells are 100% urban.'.format(
                              ncells=np.isclose(f_urban, 1.0).sum())))
        (landusef, f_urban) = deal_with_100pct_urban(landusef, f_urban)
    lu_axis = 1  # axes of landusef are [time, landuse, x, y]
    nPFT = landusef.shape[lu_axis]
    # remove all urban land
    landusef[:, (nc.ISURBAN - 1), ...] = 0.0
    # calculate the constant factor to multiply into non-water,
    # non-urban land use fractions to fill in the vacated urban areas
    # without changing the relative proportions of non-urban land
    # covers.
    idx_not_urban = np.setdiff1d(np.arange(nPFT), [nc.ISURBAN - 1])
    sum_non_urban = np.sum(landusef[:, idx_not_urban, ...], axis=lu_axis)
    if np.isnan(sum_non_urban).any():
        raise ValueError('NaNs in land use fraction sums')
    if np.isclose(sum_non_urban, 0.0).any():
        raise ValueError('Some cells are 100% urban; dealing '
                         'with these is not implemented.')
    c_adj = 1.0 + (f_urban / sum_non_urban)
    landusef *= c_adj
    # make sure LANDUSEF sums to 1.0 for all pixels
    try:
        assert(np.allclose(np.nansum(landusef, axis=lu_axis), 1.0))
    except AssertionError as e:
        print('LANDUSEF does not sum to 1.0 for all pixels')
        return()
    nc.variables['LANDUSEF'][...] = landusef

    # now reset LU_INDEX to largest remaining land use
    lu_index = landusef.argmax(axis=lu_axis) + 1
    nc.variables['LU_INDEX'][...] = lu_index

    nc.history = ("created by metgrid_exe.  Urban land use fraction set "
                  "to 0.0 and all other non-water land covers increased "
                  "proportionally.")
    nc.close()
    print('reduced urban fraction to 0.0 in {}'.format(fname_wrf))


def km_to_yatir(lon, lat):
    """calculate kilometers from (lon, lat) to center of Yatir Forest

    This is a helper function for use_yatir_parameterization()

    Yatir is at 35.052224 E, 31.345315 N

    ARGUMENTS
    lon (array-like): array of longitude values (degrees E)
    lat (array-like): array of latitude values (degrees N)
    """
    arr_shape = lon.shape
    if not (lat.shape == arr_shape):
        raise ValueError('lon, lat must have same shape')
    d_km = np.full(lon.flatten().shape, np.nan)
    yatir_lon, yatir_lat = 35.052224, 31.345315
    for i in range(lon.flatten().size):
        d_km[i] = geopy.distance.geodesic((yatir_lat, yatir_lon),
                                          (lat.flatten()[i], lon.flatten()[i])).km
    return(d_km.reshape(arr_shape))


def use_yatir_parameterization(fname_wrf, dist_cutoff=16):
    """change the landuse category for Yatir forest to 20

    changes LANDUSE all WRF pixels within dist_cutoff kilometers of
    Yatir Forest to 20.

    20 matches the new entry for Yatir Forest in VEGPARM.TBL and
    LANDUSE.TBL
    """
    nc = netCDF4.Dataset(fname_wrf, 'a')  # open in append mode
    # I'm hard coding to 22 so I can test this routine repeatedly
    # without incrementing NUM_LAND_CAT every time.
    # type(nc.NUM_LAND_CAT)() casts the 22 to the same integer type
    # that nc.NUM_LAND_CAT already is.  TODO: is there a better way to
    # code this?
    nc.NUM_LAND_CAT = type(nc.NUM_LAND_CAT)((22))
    if any([substr in fname_wrf for substr in ("wrfbdy",
                                               "wrflow",
                                               "wrfinput_d01")]):
        # for these files, only change NUM_LAND_CAT attribute
        nc.close()
        return()

    LU_YATIR = 20
    WRF_lon = nc.variables['XLONG_C'][...].squeeze()
    WRF_lat = nc.variables['XLAT_C'][...].squeeze()
    d_yatir = km_to_yatir(WRF_lon, WRF_lat)
    idx_yatir = np.argwhere(d_yatir < 5)
    YATIR_X = idx_yatir[:, 0]
    YATIR_Y = idx_yatir[:, 1]
    landuse = nc.variables['LU_INDEX'][...]
    landuse[:, YATIR_X, YATIR_Y] = LU_YATIR
    nc.variables['LU_INDEX'][...] = landuse
    print("modified LU_INDEX for Yatir in {fname_wrf}".format(fname_wrf=fname_wrf))

    # LANDUSEF specifies a fraction for all land use types in each
    # pixel.  Set all urban fraction in LANDUSEF to (1.0 -
    # LANDUSEF[LU_water]) and all other non-water fractions to 0.0.
    # need to subtract 1 from land use codes to translate them to
    # zero-based array indices
    landusef = nc.variables['LANDUSEF'][...]
    # set all LU fractions to 0.0 in Yatir pixel
    landusef[:, :, YATIR_X, YATIR_Y] = 0.0
    # set Yatir LU fraction to 1.0 in Yatir pixel
    landusef[:, LU_YATIR - 1, YATIR_X, YATIR_Y] = 1.0
    nc.variables['LANDUSEF'][...] = landusef
    nc.history = ("created by metgrid_exe.  Land use set "
                  "to Yatir in WRF pixels within {} km"
                  " of Yatir Forest").format(dist_cutoff)
    nc.close()


def make_redwood_range_urban(redwoods_mask, fname_wrf):
    """make redwood range urban using digital redwoods data
    """
    nc = netCDF4.Dataset(fname_wrf, 'a')  # open in append mode
    # set LU_INDEX to urban for all non-water pixels.  This approach
    #
    landuse = nc.variables['LU_INDEX'][...]
    landuse[np.where(redwoods_mask)] = nc.ISURBAN
    nc.variables['LU_INDEX'][...] = landuse
    print("modified LU_INDEX in {fname_wrf}".format(fname_wrf=fname_wrf))

    # LANDUSEF specifies a fraction for all land use types in each
    # pixel.  Set all urban fraction in LANDUSEF to (1.0 -
    # LANDUSEF[LU_water]) and all other non-water fractions to 0.0.
    # need to subtract 1 from land use codes to translate them to
    # zero-based array indices
    landusef = nc.variables['LANDUSEF'][...]
    redwoods_mask = np.broadcast_to(redwoods_mask,
                                    landusef.shape)
    plant_mask = np.copy(redwoods_mask)
    plant_mask[:, nc.ISURBAN - 1, ...] = False
    plant_mask[:, nc.ISWATER - 1, ...] = False
    landusef[plant_mask] = 0.0

    non_water_fraction = (
        1.0 - landusef[:, (nc.ISWATER - 1), ...])
    non_water_mask = redwoods_mask[0, (nc.ISURBAN - 1), ...]
    landusef[:, (nc.ISURBAN - 1), ...][
        non_water_mask[np.newaxis, ...]] = non_water_fraction[
            non_water_mask[np.newaxis, ...]]
    nc.variables['LANDUSEF'][...] = landusef
    nc.history = ("created by metgrid_exe.  Land use set "
                  "to urban in WRF pixels containing redwoods "
                  "according to Little (1971).  Atlas of United States "
                  "trees. Volume 1. Conifers and important hardwoods.")
    nc.close()


def make_redwood_range_urban_quick_dirty(fname_wrf):
    """'quick and dirty' method to make approximate redwood range urban
    """
    nc = netCDF4.Dataset(fname_wrf, 'a')  # open in append mode
    # set LU_INDEX to urban for all non-water pixels.  This approach
    #
    crescent_city_lat = 41.7558  # Crescent City, CA latitude
    big_sur_latitude = 35.84
    landuse = nc.variables['LU_INDEX'][...]

    # met_em* uses CLAT, wrfinput* uses XLAT.
    try:
        lat = nc.variables['XLAT'][...]
    except KeyError:
        lat = nc.variables['CLAT'][...]
    near_coast = np.zeros(lat.shape, dtype='bool')
    for i in range(landuse.shape[1]):
        land_count = 0
        for j in range(landuse.shape[2]):
            if (landuse[0, i, j] != nc.ISWATER):
                if (land_count < 4):
                    near_coast[0, i, j] = True
                land_count = land_count + 1
    redwoods_mask = ((lat <= crescent_city_lat) &
                     (lat >= big_sur_latitude) &
                     near_coast)
    landuse[np.where(redwoods_mask)] = nc.ISURBAN
    nc.variables['LU_INDEX'][...] = landuse
    print("modified LU_INDEX in {fname_wrf}".format(fname_wrf=fname_wrf))

    # LANDUSEF specifies a fraction for all land use types in each
    # pixel.  Set all urban fraction in LANDUSEF to (1.0 -
    # LANDUSEF[LU_water]) and all other non-water fractions to 0.0.
    # need to subtract 1 from land use codes to translate them to
    # zero-based array indices
    landusef = nc.variables['LANDUSEF'][...]
    redwoods_mask = np.broadcast_to(redwoods_mask,
                                    landusef.shape)
    plant_mask = np.copy(redwoods_mask)
    plant_mask[:, nc.ISURBAN - 1, ...] = False
    plant_mask[:, nc.ISWATER - 1, ...] = False
    landusef[plant_mask] = 0.0

    non_water_fraction = (
        1.0 - landusef[:, (nc.ISWATER - 1), ...])
    non_water_mask = redwoods_mask[0, (nc.ISURBAN - 1), ...]
    landusef[:, (nc.ISURBAN - 1), ...][
        non_water_mask[np.newaxis, ...]] = non_water_fraction[
            non_water_mask[np.newaxis, ...]]
    nc.variables['LANDUSEF'][...] = landusef
    nc.history = ("created by metgrid_exe.  Land use set "
                  "to urban from coastline to four pixels "
                  "inland by adjust_WRF_inputs python module.")
    nc.close()


def make_all_land_urban(fname_wrf):
    """set all land use in a real.exe wrfinput netcdf file to urban

    real.exe takes meteorological data prepared by WPS and creates
    wrfinput files.  WPS is the `WRF preprocessing system
    <http://www2.mmm.ucar.edu/wrf/users/>`.

    ARGS
    fname_wrf (character): full path to a met_em* file produced by the WPS.

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
    nc = netCDF4.Dataset(fname_wrf, 'a')  # open in append mode
    # set LU_INDEX to urban for all non-water pixels.  This approach
    #
    landuse = nc.variables['LU_INDEX'][...]
    landuse[np.where(landuse != nc.ISWATER)] = nc.ISURBAN
    nc.variables['LU_INDEX'][...] = landuse
    print("modified LU_INDEX in {fname_wrf}".format(fname_wrf=fname_wrf))
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


def reduce_soil_moisture(fname_wrf, f, soil_moist_vars=None,
                         land_sea_var='LANDMASK'):
    """multiply WPS soil moisture values in a netcdf file by a constant factor

    WPS is the `WRF preprocessing system
    <http://www2.mmm.ucar.edu/wrf/users/>`.

    ARGS
    fname_wrf (character): full path to a met_em* file produced by the WPS.
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

    nc = netCDF4.Dataset(fname_wrf, 'a')  # open in append mode
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
        print("modified {this_var} in {fname_wrf}").format(this_var=this_var,
                                                           fname_wrf=fname_wrf)
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
