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

from matplotlib.cm import get_cmap
from matplotlib.figure import Figure

from map_tools_twh.map_tools_twh import CoastalSEES_WRF_prj
from map_tools_twh.map_tools_twh import CoastalSEES_WRF_Mapper

from timutils.midpt_norm import get_discrete_midpt_cmap_norm
from timutils.colormap_nlevs import setup_colormap


class MyFig(Figure):
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

    def read_soil_layers(self, silent=False):
        """read soil layers, print to stdout
        """
        nf = netCDF4.MFDataset(self.fnames['control'])
        # ZS is soil layer midpoints
        zs = nf.variables['ZS'][0, ...]  # assume (for now) that soil
                                         # layers are time-invariant
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

    def get_layer_str(self, layer):
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
        for k, v in self.data.iteritems():
            nf = netCDF4.MFDataset(self.fnames[k])
            self.data[k] = nf[self.varname][...]
            # read latitude
            if self.lat is None:
                self.lat = nf['XLAT'][0, ...]
            else:
                if np.allclose(nf['XLAT'][0, ...],
                               self.lat,
                               equal_nan=True) is False:
                    raise RuntimeError(error_str.format(labA=self.label_A,
                                                        labB=self.label_B,
                                                        var='latitude'))
            # read longitude
            if self.lon is None:
                self.lon = nf['XLONG'][0, ...]
            else:
                if np.allclose(nf['XLONG'][0, ...],
                               self.lon,
                               equal_nan=True) is False:
                    raise RuntimeError(error_str.format(labA=self.label_A,
                                                        labB=self.label_B,
                                                        var='longitude'))
            # read land/water mask
            if self.is_land is None:
                self.is_land = nf['XLAND'][0, ...] == 1

            # read units
            if self.units is None:
                self.units = nf[self.varname].units
            elif nf[self.varname].units != self.units:
                raise RuntimeError(error_str.format(labA=self.label_A,
                                                    labB=self.label_B,
                                                    var='units'))
            # read time
            xtime = nf["XTIME"][...].astype('float')
            sim_start_time = nf.SIMULATION_START_DATE
            self.time[k] = pd.DatetimeIndex([
                datetime.datetime.strptime(sim_start_time,
                                           "%Y-%m-%d_%H:%M:%S") +
                pd.tseries.offsets.DateOffset(minutes=m) for m in xtime])

            self.longname = nf[self.varname].description
            nf.close()
        self._match_tstamps()

    def _match_tstamps(self):
        """find time corresponding time indices
        """
        idx_A = self.time[self.label_A].isin(self.time[self.label_B])
        idx_B = self.time[self.label_B].isin(self.time[self.label_A])
        self.data[self.label_A] = self.data[self.label_A][idx_A, ...]
        self.data[self.label_B] = self.data[self.label_B][idx_B, ...]
        self.time = self.time[self.label_A][idx_A]

    def mask_oceans(self):
        """mask water pixels in data
        """
        if self.is_land is not None:
            for k in self.data.keys():
                ocean_mask = np.broadcast_to(np.logical_not(self.is_land),
                                             self.data[k].shape)
                self.data[k] = ma.masked_where(ocean_mask, self.data[k])


def graphics(vd, t_idx=0, layer=None, fig_type='png'):
    """plot contours of WRF var vals, differences from two different runs

    ARGS:
    vd (var_diff object): values for one variable from two different
       model runs
    t_idx (int): time stamp index (in [0, number of time stamps])
       layer (int): the vertical layer to be plotted.
    fig_type ({"png"}|"pdf"): type of image to create
    """

    t0 = datetime.datetime.now()
    # construct index into data
    if layer is None:
        idx = np.s_[t_idx, ...]
        layer_id = ""
    else:
        idx = np.s_[t_idx, layer, ...]
        layer_id = "lay{}_".format(layer)

    fname = os.path.join('/global/homes/t/twhilton',
                         'plots', 'Summen',
                         "{varname}_{layer_id}diff_maps_{tstamp}.{ext}".format(
                             varname=vd.varname,
                             layer_id=layer_id,
                             tstamp=vd.time[t_idx].strftime('%Y-%m-%d_%H%M'),
                             ext=fig_type))
    print('plotting {}'.format(fname))

    # initialize figure, axes
    nplots = 4
    fig = MyFig(figsize=(8, 8))
    ax = [None] * nplots
    for axidx, axspec in enumerate(range(221, 225)):
        ax[axidx] = fig.add_subplot(axspec, projection=CoastalSEES_WRF_prj())
        ax[axidx].set_extent((vd.lon.min(), vd.lon.max(),
                              vd.lat.min(), vd.lat.max()))

    dmin = 0.0  # min(all_data)
    dmax = np.max(map(np.max, vd.data.values()))
    cmap, norm = setup_colormap(dmin, dmax, nlevs=10,
                                cmap=get_cmap('YlGnBu'))

    for axidx, k in enumerate(vd.data.keys()):
        print("    plot {} data - {}".format(k, str(vd.time[t_idx])))
        this_map = CoastalSEES_WRF_Mapper(ax=ax[axidx])
        # TODO: mask oceans.  Probably should have an argument to
        # control this; would probably want oceans for e.g. latent
        # heat flux, but not for soil moisture

        this_map.pcolormesh(vd.lon,
                            vd.lat,
                            vd.data[k][idx],
                            norm=norm, cmap=cmap)
        this_map.colorbar()
        this_map.ax.set_title(k)

    # plot the difference
    d = vd.data[vd.label_A][idx] - vd.data[vd.label_B][idx]
    d = ma.masked_where(np.isclose(d, 0.0), d)
    idx_max = vd.data[vd.label_A].shape[0]
    if layer is None:
        idxA = np.s_[...]
        idxB = np.s_[:idx_max, ...]
    else:
        idxA = np.s_[:, layer, ...]
        idxB = np.s_[:idx_max, layer, ...]
    d_all = vd.data[vd.label_A][idxA] - vd.data[vd.label_B][idxB]
    abs_max = np.abs((d_all.min(), d_all.max())).max()
    cmap, norm = get_discrete_midpt_cmap_norm(vmin=abs_max * -1.0,
                                              vmax=abs_max,
                                              midpoint=0.0,
                                              this_cmap=get_cmap('BrBG'))
    d_map = CoastalSEES_WRF_Mapper(ax=ax[2])
    d_map.pcolormesh(vd.lon, vd.lat, d, cmap=cmap, norm=norm)
    d_map.colorbar()
    d_map.ax.set_title("{labA} - {labB} ({units})".format(
        labA=vd.label_A,
        labB=vd.label_B,
        units=vd.units))

    d_pct = (d / vd.data[vd.label_A][idx]) * 100.0
    d_pct_all = (d_all / vd.data[vd.label_A][idxA]) * 100.0
    abs_max = np.abs((d_pct_all.min(), d_pct_all.max())).max()
    cmap, norm = get_discrete_midpt_cmap_norm(vmin=abs_max * -1.0,
                                              vmax=abs_max,
                                              midpoint=0.0,
                                              this_cmap=get_cmap('BrBG'))
    pct_map = CoastalSEES_WRF_Mapper(ax=ax[3])
    pct_map.pcolormesh(vd.lon, vd.lat, d_pct, cmap=cmap, norm=norm)
    pct_map.colorbar()
    pct_map.ax.set_title("{labA} to {labB} pct change".format(
        labA=vd.label_A,
        labB=vd.label_B))

    # Draw a title before we draw plots
    title = "{vname}, {layerid}{tstamp} UTC ({units})".format(
        vname=vd.longname,
        layerid="{}, ".format(vd.get_layer_str(layer)),
        tstamp=vd.time[t_idx].strftime('%d %b %Y %H:%M'),
        units=vd.units)
    fig.suptitle(title)
    fig.savefig(fname=fname)
    print("done ({})".format(str(datetime.datetime.now() - t0)))
    return(None)


if __name__ == "__main__":
    cscratch = os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton')
    ctl_dir = os.path.join(cscratch, 'WRFv3.9_Sensitivity',
                           'WRFv3.9_Sensitivity_Ctl', 'WRFV3',
                           'run', 'summen_sensitivity_ctl')
    dry_dir = os.path.join(cscratch, 'WRFv3.9_Sensitivity',
                           'WRFv3.9_Sensitivity_DrySoil', 'WRFV3',
                           'run', 'summen_sensitivity_drysoil')
    vd = var_diff(os.path.join(ctl_dir, 'wrfsees_ccs3pb1_ls2_d02_2009-06*'),
                  os.path.join(dry_dir, 'wrfsees_ccs3pb1_ls2_d02_2009-06*'),
                  label_A='control',
                  label_B='dry',
                  varname='SMOIS')
    read_data = True
    if read_data:
        vd.read_files()
        vd.mask_oceans()
    for this_t in range(2, 100):  # 255
        fig = graphics(vd, t_idx=this_t, layer=0)
