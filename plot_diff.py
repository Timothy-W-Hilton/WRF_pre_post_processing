"""plot contours of WRF variable from two different runs and their difference

Produces a three-panel figure showing contours overlaied on a map for
a specified WRF variable.  The three panels show the values from WRF
run A, run B, and (run A - run B).

Timothy W. Hilton, UC Merced, thilton@ucmerced.edu
"""

import numpy as np
import pandas as pd
import numpy.ma as ma
import datetime
import os
import netCDF4
import glob
from wrf import getvar, extract_times, to_np, ALL_TIMES

from matplotlib.cm import get_cmap
from matplotlib.figure import Figure
from matplotlib import colors


from map_tools_twh.map_tools_twh import CoastalSEES_WRF_prj
from map_tools_twh.map_tools_twh import CoastalSEES_WRF_Mapper

from timutils.midpt_norm import get_discrete_midpt_cmap_norm
from timutils.colormap_nlevs import setup_colormap


class MyFig(Figure):
    """Subclass of matplotlib.figure.Figure; provides one-line saving

    matplotlib.figure.Figure requires some boilerplate to save a
    figure outside of matplotlib.pyplot.  Using pyplot often doesn't
    play well with detaching and reattaching screen sessions because
    screen loses its connection to $DISPLAY.  This class allows
    one-line figure saving independent of platform and $DISPLAY.

    """
    def savefig(self, dpi=150, fname="figure.pdf"):
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
        # The size * the dpi gives the final image sys.getsizeof()
        #   a4"x4" image * 80 dpi ==> 320x320 pixel image
        canvas.print_figure(fname, dpi=dpi)


class wrf_var(object):
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

    def read_files(self, mask_land=False, mask_water=False):
        """read variable from run output
        """
        nclist = [netCDF4.Dataset(f, mode="r") for f in self.fnames]
        self.data = getvar(nclist, varname=self.varname, timeidx=ALL_TIMES)
        if self.units is None:
            self.units = self.data.units
        if self.lat is None:
            self.lat = self.data.coords['XLAT'].values
        if self.lon is None:
            self.lon = self.data.coords['XLONG'].values
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


class var_diff(object):
    def __init__(self, fname_A, fname_B, label_A, label_B, varname):
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

    def read_soil_layers(self, silent=False):
        """read soil layers, print to stdout
        """
        nf = netCDF4.MFDataset(self.fnames[self.label_A])
        # ZS is soil layer midpoints
        # assume (for now) that soil layers are time-invariant
        zs = nf.variables['ZS'][0, ...]
        # DZS is soil layer thickness
        dzs = nf.variables['DZS'][0, ...]
        depth_top = np.zeros(len(zs))
        depth_bot = np.zeros(len(zs))
        for this_lay in range(len(zs)):
            depth_top[this_lay] = zs[this_lay] - (dzs[this_lay] / 2.0)
            depth_bot[this_lay] = zs[this_lay] + (dzs[this_lay] / 2.0)
        if silent is False:
            print("soil layer {}: {:0.1f} - {:0.1f} m".format(
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
        if self.varname.upper() == "LCL":
            lcl_flag = True
            # wrf.getvar() calculates LCL as part of reading analysis
            # variable 'cape_2d'.  for cape_2d wrf.getvar() returns an
            # array of shape (4, ...) of which axis 2 is LCL.
            self.varname = 'cape_2d'
        else:
            lcl_flag = False
        for k, v in self.data.items():
            nf = netCDF4.MFDataset(self.fnames[k])
            # TODO: decide whether to keep or get rid of nf.
            wv = wrf_var(sorted(glob.glob(self.fnames[k])),
                         label=self.label_A,
                         varname=self.varname,
                         is_atm=False)
            wv.read_files()
            self.data[k] = to_np(wv.data)

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
            # read land/water mask
            if self.is_land is None:
                self.is_land = np.isclose(to_np(getvar(nf, 'XLAND')), 1.0)

            # read units
            if self.units is None:
                self.units = wv.units
            elif wv.units != self.units:
                raise RuntimeError(error_str.format(labA=self.label_A,
                                                    labB=self.label_B,
                                                    var='units'))
            # read variable description to longname
            self.longname = wv.data.description
            # read time
            xtime = extract_times(nf, ALL_TIMES)
            self.time[k] = pd.DatetimeIndex(xtime)
            # read model heights
            try:
                self.z = getvar(nf, 'z')
            except ValueError as e:
                print('unable to read Z from input file: ' + str(e))
            nf.close()
        self._match_tstamps()
        if lcl_flag:
            # extract LCL from wrf.getvar('cape_2d') output
            for k, v in self.data.items():
                LCL_IDX = self.longname.lower().split(' ; ').index('lcl')
                self.data[k] = self.data[k][LCL_IDX, ...]
            self.longname = 'LCL'
            self.units = self.units.split(';')[LCL_IDX].strip()
        print('done var_diff.read_files()')

    def _match_tstamps(self):
        """find time corresponding time indices
        """
        idx_A = self.time[self.label_A].isin(self.time[self.label_B])
        idx_B = self.time[self.label_B].isin(self.time[self.label_A])
        self.data[self.label_A] = np.take(self.data[self.label_A],
                                          np.flatnonzero(idx_A),
                                          self.var_axes['Time'])
        self.data[self.label_B] = np.take(self.data[self.label_B],
                                          np.flatnonzero(idx_B),
                                          self.var_axes['Time'])
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
        if self.data[self.label_A].ndim is 3:
            idx = np.s_[t_idx, ...]
        elif self.data[self.label_A].ndim is 4:
            idx = np.s_[t_idx, layer, ...]
        elif self.data[self.label_A].ndim is 5:
            idx = np.s_[:, t_idx, layer, ...]
        else:
            raise IndexError("data have unexpected shape")
        return(idx)

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
        d_all = (self.data[self.label_A][idxA] -
                 self.data[self.label_B][idxB])
        self.abs_max = np.abs((d_all.min(), d_all.max())).max()
        self.d_pct = (self.d / self.data[self.label_A][idx]) * 100.0
        # d_pct_all = (d_all / self.data[self.label_A][idxA]) * 100.0


def wrf_var_find_axes(wv):
    """locate vertical, horizontal, and time axes in a WRF variable

    ARGS:
    wv (xarray): XArray containing the WRF output, as read by
       wrf.getvar() or wrf.extract_vars().

    RETURNS
    a dict containing the integer values of the time, vertical,
       north-south, and east-west axes of the variable.
    """
    axes = {}
    axes.update({"Time": wv.data.dims.index('Time')})
    try:
        axes.update({"Lay": wv.data.dims.index('bottom_top')})
    except ValueError:
        print("variable has no atmospheric vertical axis")
        axes.update({"Lay": None})
        try:
            axes.update({"Lay": wv.data.dims.index('soil_layers_stag')})
        except ValueError:
            print("variable has no soil axis")
            axes.update({"Lay": None})
    axes.update({"Lon": wv.data.dims.index('west_east')})
    axes.update({"Lat": wv.data.dims.index('south_north')})
    return(axes)


class VarDiffPlotter(object):
    """plot contours of WRF var vals, differences from two different runs
    """

    def __init__(self, vd, t_idx=0, layer=None, fig_type='png',
                 domain=2, pfx=None):
        """
        Initialize a VarDiffPlotter with a Figure instance and four Axes

        ARGS:
        vd (var_diff instance): values for one variable from two
           different model runs
        t_idx (int): time stamp index (in [0, number of time stamps])
           layer (int): the vertical layer to be plotted.
        fig_type ({"png"}|"pdf"): type of image to create
        """
        self.vd = vd
        self.t_idx = t_idx
        self.layer = layer
        self.fig_type = fig_type
        self.domain = domain
        self.pfx = pfx

    def get_filename(self):
        """return string containing filename for saving plot
        """
        if self.pfx is not None:
            self.pfx = self.pfx + '_'
        self.fname = os.path.join(
            '/global/homes/t/twhilton',
            'plots', 'Summen',
            ("{pfx}{varname}_{layer_id}"
             "d{domain:02d}_diff_maps_{tstamp}.{ext}").format(
                 pfx=self.pfx,
                 varname=self.vd.varname,
                 layer_id=self._get_layer_id(),
                 domain=self.domain,
                 tstamp=self.vd.time[self.t_idx].strftime(
                     '%Y-%m-%d_%H%M'),
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
                             projection=CoastalSEES_WRF_prj())
        ax.set_extent((self.vd.lon.min(), self.vd.lon.max(),
                       self.vd.lat.min(), self.vd.lat.max()))
        this_map = CoastalSEES_WRF_Mapper(ax=ax, domain=self.domain)
        hgt = ma.masked_greater(self.vd.z[layer, ...], 400)
        this_map.pcolormesh(self.vd.lon,
                            self.vd.lat,
                            hgt,
                            vmin=0, vmax=400)
        ax.set_title('WRF height, layer {:02d} (m)'.format(layer))
        this_map.colorbar()
        fig.savefig(fname='height_lay{:02d}.png'.format(layer))

    def plot(self, cb_orientation='vertical', vmin=None, vmax=None):
        """plot contour plots for both variables, diff, pct diff

        ARGS
        cb_orientation ({"vertical"} | "horizontal"): orientation for
            the colorbars
        """
        t0 = datetime.datetime.now()

        idx = self.vd.get_tstep_idx(self.t_idx, self.layer)
        self.get_filename()
        print('plotting {}'.format(self.fname))

        # initialize figure, axes
        nplots = 4
        fig = MyFig(figsize=(8, 8))
        ax = [None] * nplots
        for axidx, axspec in enumerate(range(221, 225)):
            ax[axidx] = fig.add_subplot(axspec,
                                        projection=CoastalSEES_WRF_prj())
            ax[axidx].set_extent((self.vd.lon.min(), self.vd.lon.max(),
                                  self.vd.lat.min(), self.vd.lat.max()))
        if vmin is None:
            dmin = 0.0  # min(all_data)
        else:
            dmin = vmin
        if vmax is None:
            dmax = np.max(list(map(np.max, self.vd.data.values())))
        else:
            dmax = vmax
        cmap = colors.ListedColormap(['#a6cee3', 'grey'])
        bounds = [-0.0001, 0.5, 1.0001]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        for axidx, k in enumerate(self.vd.data.keys()):
            print("    plot {} data - {}".format(
                k, str(self.vd.time[self.t_idx])))
            this_map = CoastalSEES_WRF_Mapper(ax=ax[axidx], domain=self.domain)

            this_map.pcolormesh(self.vd.lon,
                                self.vd.lat,
                                self.vd.data[k][idx],
                                norm=norm, cmap=cmap)
            this_map.colorbar(orientation=cb_orientation)
            this_map.cb.set_ticks([0.25, 0.75])
            this_map.cb.set_ticklabels(['no fog', 'fog'])
            if cb_orientation is "horizontal":
                this_map.cb.ax.set_xticklabels(
                    this_map.cb.ax.get_xticklabels(),
                    rotation=-60)
            this_map.ax.set_title(k)

        # plot the difference
        self.vd.calc_diff(idx, self.layer)
        if vmax is None:
            vmax = self.vd.abs_max
        # the colorbar should always be symmetict about 0.  Make the
        # choice to use vmax to define the range if both vmin and vmax
        # are specified.  Specifying both vmin and vmax is useful for
        # e.g. LCL, where the two LCL maps should go 0 to some max
        # value, but the *difference* plots should still be symmetric
        # about 0.0.
        vmin = vmax * -1.0
        cmap = colors.ListedColormap(['#fc8d62', '#8da0cb'])
        bounds = [-1.0001, 0.0, 1.0001]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        d_map = CoastalSEES_WRF_Mapper(ax=ax[2], domain=self.domain)
        d_map.pcolormesh(self.vd.lon, self.vd.lat, self.vd.d,
                         cmap=cmap, norm=norm)
        d_map.colorbar(orientation=cb_orientation)
        d_map.cb.set_ticks([-0.5, 0.5])
        d_map.cb.set_ticklabels(['ctl no fog,\nurb. fog',
                                 'ctl fog,\nurb. no fog'])
        if cb_orientation is "horizontal":
            d_map.cb.ax.set_xticklabels(
                d_map.cb.ax.get_xticklabels(),
                rotation=-60)
        d_map.ax.set_title("{labA} - {labB} difference".format(
            labA=self.vd.label_A,
            labB=self.vd.label_B))

        # plot the pct difference
        abs_max = 150  # np.abs((d_pct_all.min(), d_pct_all.max())).max()
        cmap, norm = get_discrete_midpt_cmap_norm(vmin=abs_max * -1.0,
                                                  vmax=abs_max,
                                                  midpoint=0.0,
                                                  this_cmap=get_cmap('cool'))
        cmap = colors.ListedColormap(['#fc8d62', '#8da0cb'])
        bounds = [-1.0001, 0.0, 1.0001]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        pct_map = CoastalSEES_WRF_Mapper(ax=ax[3],
                                         domain='bigbasin',
                                         res='10m')
        pct_map.pcolormesh(self.vd.lon, self.vd.lat,
                           self.vd.d, cmap=cmap, norm=norm)
        pct_map.colorbar(orientation=cb_orientation)
        pct_map.cb.set_ticks([-0.5, 0.5])
        pct_map.cb.set_ticklabels(['ctl no fog,\nurb. fog',
                                   'ctl fog,\nurb. no fog'])
        if cb_orientation is "horizontal":
            pct_map.cb.ax.set_xticklabels(
                pct_map.cb.ax.get_xticklabels(),
                rotation=-60)
        pct_map.ax.set_title("{labA} - {labB} difference".format(
            labA=self.vd.label_A,
            labB=self.vd.label_B,
            units=self.vd.units))

        # Draw a title before we draw plots
        title = "{vname}, {layerid}{tstamp} UTC ({units})".format(
            vname=self.vd.longname,
            layerid="{}, ".format(self._get_layer_id()),
            tstamp=self.vd.time[self.t_idx].strftime('%d %b %Y %H:%M'),
            units=self.vd.units)
        fig.suptitle(title)
        fig.savefig(fname=self.fname)
        print("done ({})".format(str(datetime.datetime.now() - t0)))
        return(None)
