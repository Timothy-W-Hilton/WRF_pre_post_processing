"""plot contours of WRF variable from two different runs and their difference

Produces a three-panel figure showing contours overlaied on a map for
a specified WRF variable.  The three panels show the values from WRF
run A, run B, and (run A - run B).

Timothy W. Hilton, UC Santa Cruz, twhilton@ucsc.edu
"""

# TODO: consider moving timestamp comparison earlier in
#    var_diff.read_file().  The comparison is currently last, which #
#    forces read_files to read *all* of both WRF runs' output files
#    even # if one run is much shorter than the other.  Alternatively,
#    could # pare down the files passed to var_diff in the driver
#    script...
# TODO: improve reading speed of netCDF files.  reading one variable
#    from one month of WRF files concatenated with ncrcat takes about
#    17 seconds.  ncrcat takes 1:35 to assemble the combined file.
#    UPDATE: not sure it's possible to read files faster than it
#    already is.  Current code takes about two minutes, close to (1:35
#    + 0:17).
# TODO: improve plotting speed by updating the data on the cartopy map
#    rather than replotting the map anew for every timestep.
# TODO: merge fog/no fog plotting and colorbars into this branch.
#    UPDATE: This is probably superceded by the now-implemented fog
#    percentage plot.
# TODO: move sum_layers from driver to plot_diff.py.  Probably rename
#    to aggregate_time or something, because it does averages as well
#    as sums and now works on time axis, not vertical axis.
# TODO: maybe revise fog base height to allow it to find base > 400 m?
# TODO: add legend to diff panels to indicate that white ~= 0.0?

import numpy as np
import pandas as pd
import numpy.ma as ma
from scipy.stats import norm
import datetime
import os
import netCDF4
import re
import socket
from xarray import DataArray
from wrf import getvar, extract_times, to_np, ALL_TIMES
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader  # for redwoods range

from matplotlib.cm import get_cmap
from matplotlib.figure import Figure

from map_tools_twh.map_tools_twh import CoastalSEES_WRF_prj
from map_tools_twh.map_tools_twh import CoastalSEES_WRF_Mapper
from map_tools_twh.map_tools_twh import get_IGBP_modMODIS_21Category_PFTs_table

from timutils.midpt_norm import get_discrete_midpt_cmap_norm
from timutils.colormap_nlevs import setup_colormap
from timutils.std_error import calc_neff


def test_ax_min():
    """test that ax_min is functioning properly
    """
    arr = np.array([[[True,  True],
                     [False,  True]],
                    [[True, False],
                     [False,  True]],
                    [[False,  True],
                     [False, False]],
                    [[False, False],
                     [False,  True]]], dtype=bool)
    correct = np.array([[0,  0],
                        [-1,  0]])
    assert(np.array_equal(ax_max(arr, 0), correct))


def ax_max(arr, axis):
    """find the index of the maximum value along an axis of an array
    """
    idx = arr.argmax(axis=axis)
    idx.data[np.all(arr.data == False, axis=axis)] = -1
    return(idx)


class MyFig(Figure):
    """Subclass of matplotlib.figure.Figure; provides one-line saving

    matplotlib.figure.Figure requires some boilerplate to save a
    figure outside of matplotlib.pyplot.  Using pyplot often doesn't
    play well with detaching and reattaching screen sessions because
    screen loses its connection to $DISPLAY.  This class allows
    one-line figure saving independent of platform and $DISPLAY.

    """
    def savefig(self, dpi=150, fname="figure.pdf", **kwargs):
        """save the object's map to an image file
        """
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.backends.backend_pdf import FigureCanvasPdf
        if fname.endswith((".pdf", ".PDF")):
            canvas = FigureCanvasPdf(self)
        elif fname.endswith((".png", ".PNG")):
            canvas = FigureCanvasAgg(self)
        else:
            raise(IOError(('unrecognized figure type.'
                           '  pdf or png are supported')))
        self.set_figwidth(27.0)
        self.set_figheight(7.5)
        # The size * the dpi gives the final image sys.getsizeof()
        #   a4"x4" image * 80 dpi ==> 320x320 pixel image
        canvas.print_figure(fname, dpi=dpi, **kwargs)


class wrf_var(object):
    """class to parse WRF output files
    """
    def __init__(self, fnames, label, varname, is_atm):
        self.fnames = fnames
        self.label = label
        self.varname = varname
        self.is_atm = None
        self.longname = None
        self.units = None
        self.time = None
        self.lat = None
        self.lon = None
        self.data = None
        self.is_land = None
        self.z = None  # height above sea level (meters)
        self.lcl_flag = False
        self.fog_present_flag = False
        self.fog_pct_flag = False
        self.fog_base_height_flag = False

    def read_land_water_mask(self):
        """read land/water mask from WRF netcdf files

        RETURNS:
        array containing True for land pixels, False for water pixels
        """
        nclist = [netCDF4.Dataset(f, mode="r") for f in self.fnames]
        xland = getvar(nclist, 'XLAND', timeidx=ALL_TIMES, meta=False)
        land_value = 1.0  # land pixels are set to 1.0, water to 2.0
        return(np.isclose(xland, land_value))

    def read_soil_layers(self, silent=False):
        """read soil layers, optionally print to stdout
        """
        nclist = [netCDF4.Dataset(f, mode="r") for f in self.fnames]
        # ZS is soil layer midpoints
        zs = getvar(nclist, 'ZS', timeidx=ALL_TIMES, meta=False)
        # DZS is soil layer thickness
        dzs = getvar(nclist, 'DZS', timeidx=ALL_TIMES, meta=False)
        depth_top = zs - (dzs / 2.0)
        depth_bot = zs + (dzs / 2.0)
        if silent is False:
            for this_lay in zs.shape[1]:
                print("soil layer {}: {:0.1f} - {:0.1f} m".format(
                    this_lay, depth_top[this_lay], depth_bot[this_lay]))
        for this_nc in nclist:
            this_nc.close()
        return({'top': depth_top, 'bot': depth_bot})

    def get_atm_layer_str(self, layer):
        """return a string describing atmosphere height above sea level
        """
        return("{} m ASL".format(self.z[layer, ...].mean()))

    def get_soil_layer_str(self, layer, t_idx=0):
        """return a string describing soil layer depth

        returns a string the format "T - B m", with T the depth of the top
        of the layer and B the depth at the bottom.

        ARGS:
        layer (int): index of the soil layer in the WRF netCDF data
        """
        if layer is None:
            return("")
        else:
            layer_depths = self.read_soil_layers(silent=True)
            return("{:0.1f} - {:0.1f} m".format(
                layer_depths['top'][t_idx, layer],
                layer_depths['bot'][t_idx, layer]))

    def get_fog_base_height(self, z_threshold=400, q_threshold=0.05):
        """find lowest vertical level in every column with qc >= 0.05 g / kg

        Implements the definition of fog from O'Brien et al 2013: fog
        is defintied is present in a horizontal gridcell if any layer
        at or below 400 m has liquid water content (qc) >= 0.05 g /
        kg.

        ARGS:
        z_threshold (int): meters above sea level to test for fog.
           Default is 400 (from O'Brien et al. (2013).
        q_threshold (int): cloud liquid water content (in g H2O / kg
           dry air) at which to consider the air "foggy". Default is
           0.05 (from O'Brien et al. (2013).

        REFERENCES

        O'Brien, T. A., L. C. Sloan, P. Y. Chuang, I. C. Faloona, and
        J. A. Johnstone (2013), Multidecadal simulation of coastal fog
        with a regional climate model, Climate Dynamics, 40(11-12),
        2801-2812, doi:10.1007/s00382-012-1486-x.
        """
        self.is_foggy_obrien_2013_3D(z_threshold, q_threshold)
        vertical_axis = 1  # axes are (0=time, 1=vertical, 2=x, 3=y)
        # self.data = ma.masked_less(ax_max(self.data, axis=vertical_axis), 0)
        zidx = ax_max(self.data, axis=vertical_axis)
        fogbase_height = DataArray(data=np.empty(zidx.shape,
                                                 dtype=float),
                                   coords=zidx.coords,
                                   dims=zidx.dims,
                                   name='fogbase_height')
        fogbase_height.data[:] = np.nan
        j, k = np.meshgrid(np.arange(self.z.shape[1]),
                           np.arange(self.z.shape[2]))
        fogbase_height.data[:, j, k] = self.z.data[zidx.data[:, j, k], j, k]
        fogbase_height.data[zidx < 0] = np.nan
        self.data = fogbase_height
        self.longname = 'fog base height'
        self.units = 'm'

    def is_foggy_obrien_2013_3D(self, z_threshold=400, q_threshold=0.05):
        """find near-surface grid cells with qc >= 0.05 g / kg

        Implements the definition of fog from O'Brien et al 2013: fog
        is defintied is present in a horizontal gridcell if any layer
        at or below 400 m has liquid water content (qc) >= 0.05 g /
        kg.

        ARGS:
        z_threshold (int): meters above sea level to test for fog.
           Default is 400 (from O'Brien et al. (2013).
        q_threshold (int): cloud liquid water content (in g H2O / kg
           dry air) at which to consider the air "foggy". Default is
           0.05 (from O'Brien et al. (2013).

        REFERENCES

        O'Brien, T. A., L. C. Sloan, P. Y. Chuang, I. C. Faloona, and
        J. A. Johnstone (2013), Multidecadal simulation of coastal fog
        with a regional climate model, Climate Dynamics, 40(11-12),
        2801-2812, doi:10.1007/s00382-012-1486-x.
        """
        t0 = datetime.datetime.now()
        print('is_foggy_obrien_2013_3D()', end='')

        gram_per_kg = 1e-3
        q_threshold = q_threshold * gram_per_kg
        cell_is_foggy = np.zeros(self.data.shape, dtype='bool')
        cell_is_foggy[(self.data >= q_threshold) &
                      (self.z <= z_threshold)] = True
        self.data.data = cell_is_foggy
        self.longname = 'fog_present_3D'
        self.units = 'Boolean'
        print(' done is_foggy_obrien_2013_3D ({})'.format(
            datetime.datetime.now() - t0))

    def is_foggy_obrien_2013_2D(self, z_threshold=400, q_threshold=0.05):
        """find horizontal cells where any layer <= 400 m has qc >= 0.05 g / kg

        Implements the definition of fog from O'Brien et al 2013: fog
        is defintied is present in a horizontal gridcell if any layer
        at or below 400 m has liquid water content (qc) >= 0.05 g /
        kg.

        ARGS:
        z_threshold (int): meters above sea level to test for fog.
           Default is 400 (from O'Brien et al. (2013).
        q_threshold (int): cloud liquid water content (in g H2O / kg
           dry air) at which to consider the air "foggy". Default is
           0.05 (from O'Brien et al. (2013).

        REFERENCES

        O'Brien, T. A., L. C. Sloan, P. Y. Chuang, I. C. Faloona, and
        J. A. Johnstone (2013), Multidecadal simulation of coastal fog
        with a regional climate model, Climate Dynamics, 40(11-12),
        2801-2812, doi:10.1007/s00382-012-1486-x.

        """
        t0 = datetime.datetime.now()
        print('is_foggy_obrien_2013_2D()', end='')
        self.longname = 'fog_present_2D'
        self.units = 'Boolean'
        self.is_foggy_obrien_2013_3D(z_threshold, q_threshold)
        vertical_axis = 1  # axes are (0=time, 1=vertical, 2=x, 3=y)
        self.data = self.data.any(axis=vertical_axis)
        print(' done is_foggy_obrien_2013_2D ({})'.format(
            datetime.datetime.now() - t0))

    def read_files(self, mask_land=False, mask_water=False):
        """read variable from run output
        """
        t0 = datetime.datetime.now()
        if self.varname.upper() == "LCL":
            self.lcl_flag = True
            # wrf.getvar() calculates LCL as part of reading analysis
            # variable 'cape_2d'.  for cape_2d wrf.getvar() returns an
            # array of shape (4, ...) of which axis 2 is LCL.
            self.varname = 'cape_2d'
        elif self.varname.lower() == 'fogpresent':
            self.fog_present_flag = True
            # presence of fog is calculated from QCLOUD
            self.varname = 'QCLOUD'
        elif self.varname.lower() == 'fogpct':
            self.fog_pct_flag = True
            # presence of fog is calculated from QCLOUD
            self.varname = 'QCLOUD'
        elif self.varname.lower() == 'fogbase':
            self.fog_base_height_flag = True
            self.varname = 'QCLOUD'
        print('start netCDF4.Dataset()', end='')
        t0_Dataset = datetime.datetime.now()
        nclist = [netCDF4.Dataset(f, mode="r") for f in self.fnames]
        print('done netCDF4.Dataset() ({})'.format(
            datetime.datetime.now() - t0_Dataset))
        print('start getvar()', end='')
        t0_getvar = datetime.datetime.now()
        self.data = getvar(nclist,
                           varname=self.varname,
                           timeidx=ALL_TIMES,
                           squeeze=False)
        print('done getvar() ({})'.format(
            datetime.datetime.now() - t0_getvar))
        if self.units is None:
            self.units = self.data.units
        if self.lat is None:
            try:
                self.lat = self.data.coords['XLAT'].values
            except KeyError as e:
                print('XLAT not present, using XLAT_M')
                self.lat = self.data.coords['XLAT_M'].values
        if self.lon is None:
            try:
                self.lon = self.data.coords['XLONG'].values
            except KeyError as e:
                print('XLONG not present, using XLONG_M')
                self.lat = self.data.coords['XLONG_M'].values
        # read time
        xtime = extract_times(nclist, ALL_TIMES)
        self.time = pd.DatetimeIndex(xtime)
        self.longname = self.data.description
        self.z = getvar(nclist, 'z')
        for this_nc in nclist:
            this_nc.close()
        if mask_land or mask_water:
            m = self.read_land_water_mask()
            if self.data.ndim == 4:
                m = m[:, np.newaxis, ...]
            m = np.broadcast_to(m, self.data.values.shape)
            if mask_water:
                m = np.logical_not(m)
            self.data.values = ma.masked_where(self.data.values, m)
        if self.lcl_flag:
            # extract LCL from wrf.getvar('cape_2d') output
            LCL_IDX = self.longname.lower().split(' ; ').index('lcl')
            self.varname = 'LCL'
            self.data = self.data[LCL_IDX, ...]
            self.longname = 'lifting condensation level'
            self.units = self.units.split(';')[LCL_IDX].strip()
        elif self.fog_present_flag:
            self.is_foggy_obrien_2013_2D()
        elif self.fog_pct_flag:
            self.get_fog_pct()
        elif self.fog_base_height_flag:
            self.get_fog_base_height()
        print('done wrf_var.read_files() ({})'.format(
            datetime.datetime.now() - t0))

    def get_fog_pct(self, z_threshold=400, q_threshold=0.05):
        """find proportion of time each horizontal cell contains fog

        Implements the definition of fog from O'Brien et al 2013: fog
        is defintied is present in a horizontal gridcell if any layer
        at or below 400 m has liquid water content (qc) >= 0.05 g /
        kg.

        ARGS:
        z_threshold (int): meters above sea level to test for fog.
           Default is 400 (from O'Brien et al. (2013).
        q_threshold (int): cloud liquid water content (in g H2O / kg
           dry air) at which to consider the air "foggy". Default is
           0.05 (from O'Brien et al. (2013).

        REFERENCES

        O'Brien, T. A., L. C. Sloan, P. Y. Chuang, I. C. Faloona, and
        J. A. Johnstone (2013), Multidecadal simulation of coastal fog
        with a regional climate model, Climate Dynamics, 40(11-12),
        2801-2812, doi:10.1007/s00382-012-1486-x.

        """

        t0 = datetime.datetime.now()
        print('start get_fog_pct()', end='')
        self.is_foggy_obrien_2013_2D(z_threshold, q_threshold)
        time_axis = 0  # axes are (0=time, 2=x, 3=y)
        n_tsteps = self.data.shape[time_axis]
        pct = (self.data.sum(axis=time_axis) / n_tsteps) * 100.0
        pct = pct.expand_dims("Time", 0)
        pct = pct.assign_coords(Time=self.data.coords['Time'][0].data)
        self.data = pct
        self.varname = 'fogpct'
        self.longname = 'fog frequency (time)'
        self.units = 'percent'
        print(' done get_fog_pct ({})'.format(datetime.datetime.now() - t0))


class var_diff(object):
    def __init__(self, fname_A=None, fname_B=None,
                 label_A=None, label_B=None,
                 varname=None,
                 ncfile=None):
        """class constructor: instantiates a var_diff object.

        ncfile specifies the full path to a previously saved netCDF
        file created by var_diff.to_netcdf()

        varname may be any valid varname for `wrf.getvar()`, as well
        as "LCL" (liftring condensation level") or "fogpresent".
        fogpresent calculates a boolean indicating whether a
        horizontal cell contains or does not contain fog.  A cell is
        considered foggy if any vertical level below 400 m above
        ground level has a cloud water content greater than or equal
        to 0.05 g / kg (O'Brien et al., 2013).

        REFERENCES

        O'Brien, T. A., L. C. Sloan, P. Y. Chuang, I. C. Faloona, and
        J. A. Johnstone (2013), Multidecadal simulation of coastal fog
        with a regional climate model, Climate Dynamics, 40(11-12),
        2801-2812, doi:10.1007/s00382-012-1486-x.

        """
        if ncfile is not None:
            nc = netCDF4.Dataset(ncfile, 'r')
            print('loading data from {}'.format(ncfile))
            self.varname = nc.varname
            self.label_A, self.label_B = nc.groups.keys()
            self.units = nc.units
            self.lat = nc.variables['lat'][...]
            self.lon = nc.variables['lon'][...]
            self.longname = nc.varname
            self.time = (
                pd.TimedeltaIndex(nc.variables['time'][...], unit='s') +
                pd.Timestamp('1970-01-01 00:00:00')
            )
            self.is_land = None
            self.z = None
            self.data = {
                self.label_A: nc.groups[self.label_A].variables[
                    self.varname][...],
                self.label_B: nc.groups[self.label_B].variables[
                    self.varname][...]
            }
            dim_names = nc.groups[self.label_A].variables[
                self.varname].dimensions
            self.var_axes = list(dim_names)
            self.var_axes[self.var_axes.index('x')] = 'Lon'
            self.var_axes[self.var_axes.index('y')] = 'Lat'
            self.var_axes[self.var_axes.index('time')] = "Time"
            # read difference p-value if it was written to netCDF
            if 'p' in nc.variables:
                self.p = nc.variables['p'][...]
            else:
                self.p = None
            self.insignificant_mask = None
            nc.close()
        else:
            if None in [fname_A, fname_B, label_A, label_B, varname]:
                raise TypeError(('must specify either a netCDF file or '
                                 'all of fname_A, fname_B, label_A, label_B, '
                                 'varname'))

            self.fnames = {label_A: fname_A, label_B: fname_B}
            self.label_A = label_A
            self.label_B = label_B
            self.varname = varname
            self.longname = None
            self.units = None
            self.time = {label_A: None, label_B: None}
            self.lat = None
            self.lon = None
            self.data = {label_A: None, label_B: None}
            self.is_land = None
            self.z = None  # height above sea level (meters)
            self.insignificant_mask = None
            self.p = None

    def read_soil_layers(self, silent=False):
        """read soil layers, print to stdout
        """
        first_fname = sorted(glob.glob(self.fnames[self.label_A]))[0]
        nf = netCDF4.Dataset(first_fname)
        # ZS is soil layer midpoints
        # assume (for now) that soil layers are time-invariant
        zs = nf.variables['ZS'][...].squeeze()
        # DZS is soil layer thickness
        dzs = nf.variables['DZS'][...].squeeze()
        depth_top = zs - (dzs / 2.0)
        depth_bot = zs + (dzs / 2.0)
        if silent is False:
            for this_lay in range(len(zs)):
                print("soil layer {}: {:0.3f} - {:0.3f} m".format(
                this_lay, depth_top[this_lay], depth_bot[this_lay]))
        nf.close()
        return({'top': depth_top, 'bot': depth_bot})

    def get_atm_layer_str(self, layer):
        """return a string describing atmosphere height above sea level
        """
        return("{} m ASL".format(self.z[layer, ...].mean()))

    def get_soil_layer_str(self, layer):
        """return a string describing soil layer depth

        returns a string the format "T - B m", with T the depth of the top
        of the layer and B the depth at the bottom.

        ARGS:
        layer (int): index of the soil layer in the WRF netCDF data
        """
        if layer is None:
            return("")
        else:
            layer_depths = self.read_soil_layers(silent=True)
            return("{:0.1f} - {:0.1f} m".format(layer_depths['top'][layer],
                                                layer_depths['bot'][layer]))

    def read_files(self):
        """read variable from run A output, run B output
        """
        error_str = '{labA} {var} differs from {labB} {var}'
        print('begin var_diff.read_files()')

        for k, v in self.data.items():
            t0 = datetime.datetime.now()
            print('starting to read {}'.format(self.varname))
            search_directory = os.path.dirname(self.fnames[k])
            re_pat = os.path.basename(self.fnames[k])
            # RE search on file basenames
            matching_files = [f for f in filter(re.compile(re_pat).search,
                                                os.listdir(search_directory))]
            # put the full directory path back in
            matching_files = [os.path.join(search_directory, f)
                              for f in matching_files]
            matching_files = sorted(matching_files)
            wv = wrf_var(matching_files,
                         label=self.label_A,
                         varname=self.varname,
                         is_atm=False)
            wv.read_files()
            print('done reading {} ({})'.format(self.varname,
                                                datetime.datetime.now() - t0))
            self.data[k] = to_np(wv.data)
            t0 = datetime.datetime.now()
            print('done {} to_np ({})'.format(
                k, datetime.datetime.now() - t0))
            # locate variable dimensions - they vary from variable to
            # variable.  e.g. Time is not always the same array axis.
            self.var_axes = wrf_var_find_axes(wv)
            # read latitude
            if self.lat is None:
                self.lat = wv.lat
            else:
                if np.allclose(wv.lat,
                               self.lat,
                               equal_nan=True) is False:
                    raise RuntimeError(error_str.format(labA=self.label_A,
                                                        labB=self.label_B,
                                                        var='latitude'))
            # read longitude
            if self.lon is None:
                self.lon = wv.lon
            else:
                if np.allclose(wv.lon,
                               self.lon,
                               equal_nan=True) is False:
                    raise RuntimeError(error_str.format(labA=self.label_A,
                                                        labB=self.label_B,
                                                        var='longitude'))
            # read units
            if self.units is None:
                self.units = wv.units
            elif wv.units != self.units:
                raise RuntimeError(error_str.format(labA=self.label_A,
                                                    labB=self.label_B,
                                                    var='units'))
            # read variable description to longname
            self.longname = wv.longname
            # set time
            self.time[k] = pd.DatetimeIndex(np.array(wv.data.Time.values,
                                                     ndmin=1))
            # maybe read land/sea mask from a single file instead of
            # from the multifile dataset?  It should not change
            # timestep to timestep
            # self.is_land = np.isclose(to_np(nf.variables['XLAND']), 1.0)

            # read model heights
            # try:
            #     # maybe put multidataset read here?
            #     self.z = getvar(nf, 'z')
            # except ValueError as e:
            #     print('unable to read Z from input file: ' + str(e))
            # nf.close()
        t0 = datetime.datetime.now()
        print('start _match_tstamps... ', end='')
        self._match_tstamps()
        print('done match_tstamps ({})'.format(datetime.datetime.now() - t0))
        print('done var_diff.read_files()')

    def _match_tstamps(self):
        """find time corresponding time indices
        """
        idx_A = self.time[self.label_A].isin(self.time[self.label_B])
        idx_B = self.time[self.label_B].isin(self.time[self.label_A])
        self.data[self.label_A] = np.take(self.data[self.label_A],
                                          np.flatnonzero(idx_A),
                                          self.var_axes.index('Time'))
        self.data[self.label_B] = np.take(self.data[self.label_B],
                                          np.flatnonzero(idx_B),
                                          self.var_axes.index('Time'))
        self.time = self.time[self.label_A][idx_A]

    def mask_land_or_water(self, mask_water=True):
        """mask land or water pixels in data

        mask_water ({True}|False): if true, mask water pixels.  If
           false, mask land pixels.
        """
        if self.is_land is not None:
            for k in self.data.keys():
                mask = np.broadcast_to(self.is_land,
                                       self.data[k].shape)
                if mask_water:
                    mask = np.logical_not(mask)
                self.data[k] = ma.masked_where(mask, self.data[k])

    def get_tstep_idx(self, t_idx, layer):
        """construct an index into the data to extract given time step
        """
        # construct index into data

        # initialize to ellipses for all dimensions
        idx = [slice(None)] * list(self.data.values())[0].ndim
        try:
            idx[self.var_axes.index('Lay')] = layer
        except ValueError as e:
            print("get_tstep_idx: no vertical layers in data")
        try:
            idx[self.var_axes.index('Time')] = t_idx
        except ValueError as e:
            print("get_tstep_idx: no time dimension in data")
            raise(e)
        return(tuple(idx))

    def aggregate_layers(self, vert_avg=False, extract_layer=None):
        """aggregate var_diff data for all time vertical layers

        Calculate the sum of each run's data across all vertical
        layers, or, optionally, the arithmetic mean.  Vertical layers
        could be either soil layers or atmospheric layers.

        ARGS:
        extract_layer (int): if specified, a single layer is extracted
           and vert_avg keywordis ignored.
        vert_avg (logical): if true, calculate the arithmetic mean
           across layers.  Default is False.

        """
        if 'Lay' not in self.var_axes:
            ValueError('object contains no vertical layer data')
        for k, v in self.data.items():
            if extract_layer is not None:
                print('extracting layer {}'.format(extract_layer))
                # we're keeping a single vertical layer
                idx = [slice(None)] * self.data[k].ndim
                # `slice(None)` places a ":" in layer index
                idx[self.var_axes.index('Lay')] = extract_layer
                self.data[k] = v[tuple(idx)].squeeze()
            else:
                # we're doing an average or sum of all vertical layers
                n_lays = v.shape[self.var_axes.index('Lay')]
                self.data[k] = np.nansum(v,
                                         axis=self.var_axes.index('Lay'),
                                         keepdims=False)
                if vert_avg:
                    self.data[k] = self.data[k] / n_lays
                    print('calculating vertical average')
                else:
                    print('calculating vertical sum')
        if vert_avg:   # outside loop so string is only appended once
            self.longname = self.longname + 'vertical avg'
        # adjust axis labels
        self.var_axes.pop(self.var_axes.index('Lay'))

    def aggregate_time(self, time_avg=False):
        """aggregate var_diff data for all time steps

        Calculate the sum of each run's data across all time steps,
        or, optionally, the arithmetic mean.

        ARGS:
        time_avg (logical): if true, calculate the arithmetic mean
           across time steps.  Default is False.
        """
        for k, v in self.data.items():
            # axes are [time, y, x]
            time_axis = self.var_axes.index('Time')
            n_tsteps = v.shape[time_axis]
            self.data[k] = np.nansum(v, axis=time_axis, keepdims=True)
            if time_avg:
                self.data[k] = self.data[k] / n_tsteps
        if time_avg:   # outside loop so string is only appended once
            self.longname = self.longname + ' time avg'

    def _get_p(self, adj_autocorr=True, idx=None):
        z = self.diff_means_test(adj_autocorr=adj_autocorr, idx=idx)
        vectorized_cdf = np.vectorize(lambda x: norm.cdf(x, 0.0, 1.0))
        p = vectorized_cdf(z)
        return(p)

    def get_significance_mask(self,
                              significance,
                              adj_autocorr=True,
                              idx=None):
        """get mask for statistically insignificant differences

        apply difference of means at specified level test to obtain a
        mask for statistically significant pixels in the data.

        places the mask in self.insignificant_mask

        ARGS:
        significance (float): significance level to test for.  Must be
           in range [0.0, 1/0]
        adj_autocorr (boolean): if true, adjust the effective number
           of independent samples according to Wilks 1995.  Defaut is
           True.
        """
        self.p = self._get_p(adj_autocorr, idx)
        self.insignificant_mask = self.p < significance

    def diff_means_test(self, adj_autocorr=True, idx=None):
        """run a paired difference of means test

        Run a standard paired difference of means test (e.g. Devore
        (1995) section 9.1) on the data.

        ARGS:
        adj_autocorr (boolean): if true, adjust the effective number
           of independent samples according to Wilks 1995.
        idx (array-like): indices into self.data items to include in
           the calculation.  Allows for estimating signficance for the
           differences of a subset of the data.

        REFERENCES

        Devore, J.L., 1995. Probability and Statistics for Engineering
        and the Sciences, 4th ed.  Brooks/Cole Publishing Co., Pacific
        Grove, California, USA.

        Wilks, D., 1995 Statistical Methods in the Atmospheric
        Sciences: An Introduction.  Academic Press, New York
        """
        if idx is None:
            idx = slice(None)  # if idx is not specified, use the whole array
        ax_time = 0  # time is axis 0 in the data array
        for k in self.data.keys():
            if adj_autocorr:
                # reduce the number of effectively independent data points
                # to account for temporal autocorrelation.
                n_eff = {k: calc_neff(v[idx].astype(float), dim=ax_time)
                         for k, v in self.data.items()}
            else:
                # assume all data points are independent
                n_eff = {k: v[idx].shape[ax_time]
                         for k, v in self.data.items()}
        # calculate test statistic z according to Devore (1995) section 9.1
        means = {k: np.mean(v[idx].astype(float), axis=ax_time)
                 for k, v in self.data.items()}
        vars = {k: np.var(v[idx].astype(float), axis=ax_time)
                for k, v in self.data.items()}
        numerator = means[self.label_A] - means[self.label_B]
        denominator = (np.sqrt((vars[self.label_A] / n_eff[self.label_A]) +
                               (vars[self.label_B] / n_eff[self.label_B])))
        # return infinity where denominator is zero
        z = np.divide(numerator, denominator,
                      out=np.full_like(numerator, np.inf),
                      where=np.abs(denominator) > 1e-12)
        return(z)

    def calc_diff(self, idx, layer):
        """calculate the variables' difference, pct diff, and absolute max diff

        populate fields d (difference), d_pct (percent difference) and
        abs_max (absolute maximum difference)

        TODO: implement a four-value scheme to indicate
        (both foggy) /
        (A foggy, B not foggy) /
        (A not foggy, B foggy) /
        (neither foggy)

        ARGS:
        idx (tuple of slice instances or indices, as from numpy.s_):
           index into the time step to calculate difference for

        """
        if self.data[self.label_A][idx].dtype is np.dtype('bool'):
            d = (self.data[self.label_A][idx].astype(float) -
                 self.data[self.label_B][idx].astype(float))
        else:
            d = (self.data[self.label_A][idx] -
                 self.data[self.label_B][idx])
        self.d = ma.masked_where(np.isclose(d, 0.0), d)
        idx_max = self.data[self.label_A].shape[0]
        if layer is None:
            idxA = np.s_[...]
            idxB = np.s_[:idx_max, ...]
        else:
            idxA = np.s_[:, layer, ...]
            idxB = np.s_[:idx_max, layer, ...]
        d_all = (self.data[self.label_A][idxA].astype(float) -
                 self.data[self.label_B][idxB].astype(float))
        self.abs_max = np.nanmax(np.abs((np.nanmin(d_all.data),
                                         np.nanmax(d_all.data))))
        self.d_pct = (self.d / self.data[self.label_A][idx]) * 100.0
        # d_pct_all = (d_all / self.data[self.label_A][idxA]) * 100.0

    def to_netcdf(self, fname):
        """write a netcdf file containing variables and their difference

        ARGS:
           fname (str): full path to the netcdf file to be written.
           If fname exists it will be deleted and replace.
        """
        nc = netCDF4.Dataset(fname, mode='w')
        # make a local copy of self.var_axes
        # simply calling `var_axes_ = self.var_axes` yields a view,
        # not a copy.  Then changing var_axes_ also changes
        # self.var_axes
        var_axes_ = [x for x in self.var_axes]
        for this_ax, this_dim in enumerate(self.var_axes):
            if this_dim == "Lat":
                this_dim = 'y'
                var_axes_[this_ax] = 'y'
            if this_dim == "Lon":
                this_dim = 'x'
                var_axes_[this_ax] = 'x'
            if this_dim == "Time":
                this_dim = 'time'
                var_axes_[this_ax] = 'time'
            key1 = list(self.data.keys())[0]
            nc.createDimension(this_dim,
                               self.data[key1].shape[this_ax])
        nc.createVariable('lat', np.float, ('y', 'x'))
        nc.createVariable('lon', np.float, ('y', 'x'))
        nc.createVariable('time', np.float, ('time'))
        nc.variables['lat'][...] = self.lat
        nc.variables['lon'][...] = self.lon
        # use unix convention of seconds since 1 Jan 1970 00:00 as
        # recommended by Unidata for storing time in netCDF files
        # (https://www.unidata.ucar.edu/software/netcdf/time/recs.html)
        nc.variables['time'][...] = np.array(
            (self.time - pd.Timestamp('1970-01-01 00:00:00')).total_seconds())
        nc.variables['time'].units = 'seconds since 1970-1-1'
        grpA = nc.createGroup(self.label_A)
        grpB = nc.createGroup(self.label_B)
        var_dtype = self.data[self.label_A].dtype
        if var_dtype is np.dtype('bool'):
            var_dtype = 'i1'
        grpA.createVariable(self.varname, var_dtype, var_axes_)
        grpB.createVariable(self.varname, var_dtype, var_axes_)
        grpA.variables[self.varname][...] = self.data[self.label_A][...]
        grpB.variables[self.varname][...] = self.data[self.label_B][...]
        if self.p is not None:
            nc.createVariable('p', np.float, ('y', 'x'))
            nc.variables['p'][...] = self.p
            nc.variables['p'].description = ('significance of the difference '
                                             'in variable values from a '
                                             'paired difference of means test')
        nc.varname = self.varname
        nc.units = self.units
        nc.close()

    def get_pval_timeseries(self, interval_hrs=12):
        """calculate time series of difference p-values

        rather than calculate one p-value for the difference for entire
        time series (as in var_diff.get_significance_mask), calculate
        values at specified intervals along the time series.  This can
        help determine how long a simulation is generally needed to get
        statistical significance.

        ARGS:
        interval_hrs (float): frequency of the calculated p-values
        """
        pval_times = pd.date_range(start=self.time[0],
                                   end=self.time[-1],
                                   freq='{}H'.format(interval_hrs))
        pval_times = pval_times[1:]  # skip time = 0
        idx_end = np.array(np.where([t in pval_times
                                     for t in self.time])).squeeze()
        key1 = list(self.data.keys())[0]
        pvals = np.full((idx_end.size,
                         self.data[key1].shape[1],
                         self.data[key1].shape[2]),
                        np.nan)
        mean_diff = np.full((idx_end.size,
                             self.data[key1].shape[1],
                             self.data[key1].shape[2]),
                            np.nan)

        for i, this_end in enumerate(idx_end):
            t0 = datetime.datetime.now()
            print('calculating pvals for start:{}'.format(this_end), end='...')
            pvals[i, ...] = self._get_p(adj_autocorr=True,
                                        idx=(range(this_end), slice(None)))
            self.calc_diff((range(this_end), slice(None), slice(None)), None)
            mean_diff[i, ...] = self.d.mean(axis=0)
            print(' done ({})'.format(datetime.datetime.now() - t0))
        self.pvals_series = pvals
        self.mean_diff_series = mean_diff
        self.t_pvals_series = pval_times
        self.idx_pvals_series = idx_end


def wrf_var_find_axes(wv):
    """locate vertical, horizontal, and time axes in a WRF variable

    ARGS:
    wv (xarray): XArray containing the WRF output, as read by
       wrf.getvar() or wrf.extract_vars().

    RETURNS
    a dict containing the integer values of the time, vertical,
       north-south, and east-west axes of the variable.
    """
    axes = list(wv.data.dims)
    try:
        axes[axes.index("bottom_top")] = "Lay"
    except ValueError:
        print("{} has no atmospheric vertical axis".format(wv.varname))
        try:
            axes[axes.index('soil_layers_stag')] = "Lay"
        except ValueError:
            print("{} has no soil axis".format(wv.varname))
    axes[axes.index('west_east')] = "Lon"
    axes[axes.index('south_north')] = "Lat"
    return(axes)


class VarDiffPlotter(object):
    """plot contours of WRF var vals, differences from two different runs
    """

    def __init__(self, vd, t_idx=0, layer=None, fig_type='png',
                 domain=2, pfx=None, savedir=None,
                 time_title_str=None, show_title=True,
                 map_prj=CoastalSEES_WRF_prj(),
                 mapper=CoastalSEES_WRF_Mapper):
        """
        Initialize a VarDiffPlotter with a Figure instance and four Axes

        ARGS:
        vd (var_diff instance): values for one variable from two
           different model runs
        t_idx (int): time stamp index (in [0, number of time stamps])
        layer (int): the vertical layer to be plotted.
        fig_type ({"png"}|"pdf"): type of image to create
        domain (int): WRF domain number (to be placed in filename)
        pfx (str): optional prefix to place in filename of every image
        savedir (str): optional full path to a directory in which to
           save the figure
        show_title ({True}|False): if true show plot titles and a main
           title
        """
        self.vd = vd
        self.t_idx = t_idx
        self.layer = layer
        self.fig_type = fig_type
        self.domain = domain
        self.pfx = pfx
        self.time_title_str = time_title_str
        self.show_title = show_title
        self.map_prj = map_prj
        self.mapper = mapper

        for k in self.vd.data.keys():
            if np.isnan(self.vd.data[k]).any():
                self.vd.data[k] = ma.masked_invalid(self.vd.data[k])

        if savedir is None:
            self.savedir = os.path.join('/', 'global', 'homes', 't',
                                        'twhilton', 'plots', 'Summen')
        else:
            self.savedir = savedir

    def get_filename(self):
        """return string containing filename for saving plot
        """
        if self.pfx is not None:
            if self.pfx[-1] is not '_':
                self.pfx = self.pfx + '_'

        self.fname = os.path.join(
            self.savedir,
            ("{pfx}{varname}_{layer_id}"
             "d{domain:02d}_diff_maps_{tstamp}.{ext}").format(
                 pfx=self.pfx,
                 varname=self.vd.varname,
                 layer_id=self._get_layer_id(),
                 domain=self.domain,
                 tstamp=self.time_title_str.replace(':', ''),
                 ext=self.fig_type))

    def _get_layer_id(self):
        """return string identifying the vertical layer of the plot
        """
        if self.layer is None:
            return("")
        else:
            return("lay{}_".format(self.layer))

    def contour_height(self, layer):
        """make a map of WRF height for specified layer
        """
        fig = MyFig(figsize=(8, 8))
        ax = fig.add_subplot(111,
                             projection=self.map_prj)
        ax.set_extent((self.vd.lon.min(), self.vd.lon.max(),
                       self.vd.lat.min(), self.vd.lat.max()))
        this_map = self.mapper(ax=ax, domain=self.domain)
        hgt = ma.masked_greater(self.vd.z[layer, ...], 400)
        this_map.pcolormesh(self.vd.lon,
                            self.vd.lat,
                            hgt,
                            vmin=0,
                            vmax=400,
                            edgecolors='face')
        ax.set_title('WRF height, layer {:02d} (m)'.format(layer))
        this_map.colorbar()
        fig.savefig(fname='height_lay{:02d}.png'.format(layer))

    def plot(self,
             cb_orientation='vertical',
             vmin=None,
             vmax=None,
             mask=None,
             cmap=get_cmap('YlGnBu')):
        """plot contour plots for both variables, diff, pct diff

        ARGS
        cb_orientation ({"vertical"} | "horizontal"): orientation for
            the colorbars
        vmin (float): minimum value for color scale.  If unspecified
            defaults to 0.0.
        vmax (float): maximum value for color scale.  If unspecified
            defaults to the maximum value in the difference data.
        mask (array-like): array of True/False values.  If specified,
            True values will be masked in the plot.  Useful for, for
            example, masking out statistically insignificant pixels.
        """

        # plot() is essentially a wrapper function to dispatch
        # _plot_scalar or _plot_vector
        t0 = datetime.datetime.now()
        self._plot_init(mask=mask)
        if 'u_v' in self.vd.var_axes:
            # heat maps make huge PDF files, but barbs do not
            self.fig_type = 'pdf'
            self.get_filename()
            self._plot_vector()
        else:
            self._plot_scalar(cb_orientation, vmin, vmax, mask, cmap)
        self.fig.savefig(fname=self.fname,
                         bbox_inches='tight',
                         dpi=self.fig.get_dpi())
        print('done plotting {} ({})'.format(
            self.fname,
            datetime.datetime.now() - t0))

    def _get_idx(self):
        """TODO write this docstring
        """
        return(self.vd.get_tstep_idx(self.t_idx, self.layer))

    def _plot_scalar(self,
                     cb_orientation='vertical',
                     vmin=None,
                     vmax=None,
                     mask=None,
                     cmap=get_cmap('YlGnBu')):
        """plot heatmaps for both variables, diff, pct diff

        This function does the heavy lifting for scalar variables.  It
        is dispatched by self.plot() depending on the dimensions and
        dimension labels.

        ARGS
        cb_orientation ({"vertical"} | "horizontal"): orientation for
            the colorbars
        vmin (float): minimum value for color scale.  If unspecified
            defaults to 0.0.
        vmax (float): maximum value for color scale.  If unspecified
            defaults to the maximum value in the difference data.
        mask (array-like): array of True/False values.  If specified,
            True values will be masked in the plot.  Useful for, for
            example, masking out statistically insignificant pixels.

        """
        if vmin is None:
            dmin = 0.0  # min(all_data)
        else:
            dmin = vmin
        if vmax is None:
            dmax = np.nanmax(list(map(np.nanmax, self.vd.data.values())))
        else:
            dmax = vmax
        cmap, norm = setup_colormap(
            dmin,
            dmax,
            nlevs=21,
            cmap=cmap)

        # pcolormesh expects the coordinates to describe the lower
        # left corner of each grid cell.  WRF's mass grid describes
        # the center.  Calculate lower left from the centers by
        # substractinv half each grid cell's width or height from the
        # coordinates.  This reaults in one fewer row and column in
        # the array.
        lons_ll = (self.vd.lon[:, :-1] -
                   (np.diff(self.vd.lon, axis=1) / 2))[:-1, :]
        lats_ll = (self.vd.lat[:-1, :] -
                   (np.diff(self.vd.lat, axis=0) / 2))[:, :-1]

        for axidx, k in enumerate(self.vd.data.keys()):
            print("    plot {} data - {}".format(
                k, str(self.vd.time[self.t_idx])))
            cm = self.main_maps[axidx].pcolormesh(
                lons_ll,
                lats_ll,
                self.vd.data[k][self._get_idx()][:-2, :-2],
                norm=norm,
                cmap=cmap,
                edgecolors='face')
            self.main_maps[axidx].colorbar(orientation=cb_orientation)
            if cb_orientation is "horizontal":
                self.main_maps[axidx].cb.ax.set_xticklabels(
                    self.main_maps[axidx].cb.ax.get_xticklabels(),
                    rotation=-60)
                self.main_maps[axidx].ax.set_title(k)

        # figure out colorscale bounds for difference maps
        if vmax is None:
            vmax = self.vd.abs_max
            # the colorbar should always be symmetict about 0.  Make the
            # choice to use vmax to define the range if both vmin and vmax
            # are specified.  Specifying both vmin and vmax is useful for
            # e.g. LCL, where the two LCL maps should go 0 to some max
            # value, but the *difference* plots should still be symmetric
            # about 0.0.
        vmin = vmax * -1.0
        cmap, norm = get_discrete_midpt_cmap_norm(vmin=vmin,
                                                  vmax=vmax,
                                                  midpoint=0.0,
                                                  bands_above_mdpt=10,
                                                  bands_below_mdpt=10,
                                                  this_cmap=get_cmap('RdBu'),
                                                  remove_middle_color=True)
        cm = self.d_map.pcolormesh(lons_ll,
                              lats_ll,
                              self.vd.d[:-2, :-2],
                              cmap=cmap,
                              norm=norm,
                              edgecolors='face')
        self.d_map.colorbar(orientation=cb_orientation)
        if cb_orientation is "horizontal":
            self.d_map.cb.ax.set_xticklabels(
                self.d_map.cb.ax.get_xticklabels(),
                rotation=-60)
            self.d_map.ax.set_title("{labA} - {labB} ({units})".format(
                labA=self.vd.label_A,
                labB=self.vd.label_B,
                units=self.vd.units))

        for this_map in [self.SFBay_map]:
            this_map.pcolormesh(lons_ll, lats_ll,
                                self.vd.d[:-2, :-2],
                                cmap=cmap,
                                norm=norm,
                                edgecolors='face')
            this_map.colorbar(orientation=cb_orientation)
            if cb_orientation is "horizontal":
                this_map.cb.ax.set_xticklabels(
                    this_map.cb.ax.get_xticklabels(),
                    rotation=-60)
        return(None)

    def _plot_vector(self,
                     cb_orientation='vertical',
                     vmin=None,
                     vmax=None,
                     mask=None,
                     cmap=get_cmap('YlGnBu')):
        """plot contour plots for both variables, diff, pct diff

        ARGS
        cb_orientation ({"vertical"} | "horizontal"): orientation for
            the colorbars
        vmin (float): minimum value for color scale.  If unspecified
            defaults to 0.0.
        vmax (float): maximum value for color scale.  If unspecified
            defaults to the maximum value in the difference data.
        mask (array-like): array of True/False values.  If specified,
            True values will be masked in the plot.  Useful for, for
            example, masking out statistically insignificant pixels.
        """

        barb_params = dict(zorder=20,
                           regrid_shape=20,
                           length=5,
                           sizes=dict(emptybarb=0.1, spacing=.2, height=0.5),
                           barb_increments=dict(half=0.5, full=1, flag=5))

        for axidx, k in enumerate(self.vd.data.keys()):
            print("    plot {} data - {}".format(
                k, str(self.vd.time[self.t_idx])))
            uv = self.vd.data[k][self._get_idx()]
            # needs projection specified expliclitly?
            self.main_maps[axidx].get_ax().barbs(
                self.vd.lon,
                self.vd.lat,
                uv[0, ...],
                uv[1, ...],
                transform=ccrs.PlateCarree(),
                **barb_params)
            self.main_maps[axidx].ax.set_title(k)

        # plot the difference
        self.vd.calc_diff(slice(None), self.layer)
        uv = self.vd.d[self._get_idx()]
        for this_map in [self.d_map, self.SFBay_map]:
            this_map.get_ax().barbs(self.vd.lon,
                                    self.vd.lat,
                                    uv[0, ...],
                                    uv[1, ...],
                                    transform=ccrs.PlateCarree(),
                                    **barb_params)
        this_map.get_ax().annotate(
            s='barb increments\n' + '\n'.join(
                '{} {} m/s'.format(k, v)
                for k, v in barb_params['barb_increments'].items()),
            xy=(0.5, 0.5),
            xycoords='axes fraction',
            xytext=(1.1, 0.5),
            textcoords='axes fraction',
            fontsize='xx-large')
        return(None)

    def _plot_init(self, mask=None):
        """plot contour plots for both variables, diff, pct diff

        ARGS
        mask (array-like): array of True/False values.  If specified,
            True values will be masked in the plot.  Useful for, for
            example, masking out statistically insignificant pixels.
        """
        self.get_filename()
        print('plotting {} difference maps'.format(self.vd.varname))
        # initialize figure, axes
        nplots = 4
        self.fig = MyFig(figsize=(22, 5.5))
        self.ax = [None] * nplots
        for axidx, axspec in enumerate(range(141, 141 + nplots)):
            self.ax[axidx] = self.fig.add_subplot(
                axspec, projection=self.map_prj)
            self.ax[axidx].set_extent((self.vd.lon.min(), self.vd.lon.max(),
                                       self.vd.lat.min(), self.vd.lat.max()))

        # plot the difference
        if not(hasattr(self.vd, 'd')):
            self.vd.calc_diff(self._get_idx(), self.layer)
        if mask is not None:
            self.vd.d = ma.masked_where(mask, self.vd.d)

        self.main_maps = [self.mapper(ax=self.ax[i],
                                      domain=self.domain)
                          for i in (0, 1)]
        self.d_map = self.mapper(ax=self.ax[2], domain=self.domain)

        self.SFBay_map = self.mapper(ax=self.ax[3],
                                     domain='redwoods',
                                     res='10m')
        for this_map in [self.d_map, self.SFBay_map]:
            if self.show_title:
                this_map.ax.set_title("{labA} - {labB} ({units})".format(
                    labA=self.vd.label_A,
                    labB=self.vd.label_B,
                    units=self.vd.units))
            if 'cori' in socket.gethostname():
                redwood_data_path = os.path.join(
                    '/', 'project', 'projectdirs', 'm2319', 'Data',
                    'Redwood_SequoiaSempervirens_extentNorthAmerica',
                    'data', 'commondata', 'data0')
            elif 'MacBook' in socket.gethostname():
                redwood_data_path = os.path.join(
                    '/', 'Users', 'tim', 'work', 'Data', 'Redwoods',
                    'Redwood_SequoiaSempervirens_extentNorthAmerica',
                    'data', 'commondata', 'data0')
            else:
                raise(FileNotFoundError)
            fname = os.path.join(redwood_data_path, 'sequsemp.shp')
            rw_shapes = list(shpreader.Reader(fname).geometries())
            # a number of places on the web
            # (e.g. https://github.com/SciTools/cartopy/issues/924) suggest
            # crs.epsg(3857) and crs.Mercator.GOOGLE are the same.  But the
            # redwood range plots differently.
            # crs.Mercator.GOOGLE.proj4_params
            # and crs.epsg(3857).proj4_params differ; that's another way to
            # illustrate this.
            proj = ccrs.Mercator.GOOGLE
            this_map.ax.add_geometries(geoms=rw_shapes, crs=proj,
                                       edgecolor='#d95f02', facecolor='none',
                                       alpha=0.5)

        # make a title before we draw plots
        if self.time_title_str is None:
            self.time_title_str = self.vd.time[self.t_idx].strftime(
                '%d %b %Y %H:%M UTC')
        if self.show_title:
            title = "{vname}, {tstamp} ({units})".format(
                vname=self.vd.longname,
                tstamp=self.time_title_str,
                units=self.vd.units)
            self.fig.suptitle(title)
        # self.fig.set_tight_layout(True) and self.fig.tight_layout()
        # both issue warnings (see
        # https://github.com/matplotlib/matplotlib/issues/1852,
        # https://stackoverflow.com/questions/15455029/python-matplotlib-agg-vs-interactive-plotting-and-tight-layout).
        # self.fig.set_tight_layout(True)
        # self.fig.tight_layout()
        self.fig.set_figwidth(22.0)
        self.fig.set_figheight(5.5)
        self.fig.set_dpi(300)
        return(None)
