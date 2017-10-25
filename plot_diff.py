"""plot contours of WRF variable from two different runs and their difference

Produces a three-panel figure showing contours overlaied on a map for
a specified WRF variable.  The three panels show the values from WRF
run A, run B, and (run A - run B).

Timothy W. Hilton, UC Merced, thilton@ucmerced.edu
"""

import numpy as np
import datetime
import os
import Ngl, Nio
import netCDF4

class var_diff(object):
  def __init__(self, fname_A, fname_B, label_A, label_B, varname):
    self.fnames = {label_A:fname_A, label_B:fname_B}
    self.label_A = label_A
    self.label_B = label_B
    self.varname = varname
    self.longname = None
    self.units = None
    self.time = None
    self.lat = None
    self.lon = None
    self.data = {label_A:None, label_B:None}

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
      # read units
      if self.units is None:
        self.units = nf[self.varname].units
      elif nf[self.varname].units != self.units:
        raise RuntimeError(error_str.format(labA=self.label_A,
                                            labB=self.label_B,
                                            var='units'))
      # read time
      if self.time is None:
        self.time = nf["XTIME"][...]
      # this test failes if either file contains a shorter, but
      # otherwise identical, time series.  TODO: Needs more complex
      # logic to deal with that
      # elif np.allclose(nf["XTIME"][...],
      #                  self.time,
      #                  equal_nan=True) is False:
      # raise RuntimeError(error_str.format(labA=self.label_A,
      #                                       labB=self.label_B,
      #                                       var='timestamps'))
      # read variable long name
      self.longname = nf[self.varname].description
      nf.close()

def graphics(vd, t_idx=0, layer=None):
  """plot contours of WRF variable from two different runs and their difference

  ARGS:
  vd (var_diff object): values for one variable from two different
     model runs
  t_idx (int): time stamp index (in [0, number of time stamps])
  layer (int): the vertical layer to be plotted.

  """

  # datetime.timedelta does not support type numpy.float32, so cast to
  # type float
  # TODO: get the start date from the netcdf file instead of hard-coding
  this_t = (datetime.datetime(2009, 6, 1) +
            datetime.timedelta(minutes=float(vd.time[t_idx])))

  #---Start the graphics section
  fname = os.path.join('/global/homes/t/twhilton',
                                  'plots', 'Summen',
                       "{}_diff_maps_{}".format(
                         vd.varname, this_t.strftime('%Y-%m-%d_%H%M')))
  wks_type = "png"
  wks = Ngl.open_wks(wks_type, fname)
  print('plotting {}'.format(fname))

  # Create resource list for customizing contour over maps
  res                        = Ngl.Resources()

  res.nglMaximize            = False
  res.nglFrame               = False
  res.nglDraw                = False

  res.mpProjection           = "CylindricalEquidistant"
  res.mpOutlineOn            = True

  #---Zoom in on plot
  res.mpLimitMode = "LatLon"    # Limit the map view.
  res.mpMinLatF     = vd.lat.min()
  res.mpMaxLatF     = vd.lat.max()
  res.mpMinLonF     = vd.lon.min()
  res.mpMaxLonF     = vd.lon.max()

  res.mpPerimOn              = False
  res.mpGridAndLimbOn        = False
  res.mpDataBaseVersion = "MediumRes"
  res.pmTickMarkDisplayMode  = "Never"

  res.cnFillMode = "RasterFill"
  res.cnFillOn               = True
  res.cnLinesOn              = False
  res.cnLineLabelsOn         = False
  res.cnLevelSelectionMode   = "ManualLevels"
  res.lbLabelBarOn           = False

  res.sfXArray               = vd.lon
  res.sfYArray               = vd.lat

  #
  # Loop 9 times and create 9 dummy plots; each group
  # of 3 plots has the same color map.
  #
  nplots = 3
  if layer is None:
    idx = np.s_[t_idx, ...]
  else:
    idx = np.s_[t_idx, layer, ...]
  plots  = []

  all_data = np.concatenate((vd.data[vd.label_A].flatten(),
                             vd.data[vd.label_B].flatten()))
  dmin = 0.0  # min(all_data)
  dmax = max(all_data)

  for k in vd.data.keys():
    print("    plot {} data - {}".format(k, str(this_t)))

    res.cnFillPalette = "WhiteBlue"
    res.cnLevelSelectionMode = "ManualLevels"
    res.cnMinLevelValF         = dmin
    res.cnMaxLevelValF         = dmax
    nlevels = 10
    res.cnLevelSpacingF        = (dmax - dmin) / nlevels
    res.tiMainString  = k
    # TODO: mask oceans.  Probably should have an argument to control this; would probably want oceans for e.g. latent heat flux, but not for soil moisture
    plots.append(Ngl.contour_map(wks, vd.data[k][idx], res))

  # plot the difference
  # TODO: unified scale across all time steps
  d = vd.data[vd.label_A][idx] - vd.data[vd.label_B][idx]
  idx_max = vd.data[vd.label_A].shape[0]
  d_all = vd.data[vd.label_A] - vd.data[vd.label_B][:idx_max, ...]
  abs_max = np.abs((d_all.min(), d_all.max())).max()
  nlevs = 10
  res.cnFillPalette = "BrownBlue12"
  res.cnLevelSelectionMode = "ManualLevels"
  res.cnLevelSpacingF      = (abs_max * 2) / nlevs
  res.cnMinLevelValF         = abs_max * -1
  res.cnMaxLevelValF         = abs_max
  res.tiMainString  = "{labA} - {labB} ({units})".format(labA=vd.label_A,
                                                         labB=vd.label_B,
                                                         units=vd.units)
  plots.append(Ngl.contour_map(wks, d, res))

  # Resources for panelling
  pres                  = Ngl.Resources()
  pres.nglFrame         = False
  pres.nglPanelLabelBar = True

  # Calculate start Y position for first row of plots
  height = 0.45            # we know this will be height of small plots
  extra  = 1.0-(2*height)
  top    = 1.0-(extra/2.)

  # Draw a title before we draw plots
  title = "{vname}, {tstamp} UTC ({units})".format(
    vname=vd.longname,
    tstamp=this_t.strftime('%d %b %Y %H:%M'),
    units=vd.units)
  txres               = Ngl.Resources()
  txres.txJust        = "BottomCenter"
  txres.txFontHeightF = 0.02
  Ngl.text_ndc(wks,title,0.5,top+0.01,txres)

  # Loop across plots and panel them on one page
  for n in range(0,2):
    # Define location in a unit square for each set of plots.
    pres.nglPanelTop    = top-(n*height)
    pres.nglPanelBottom = top-((n+1)*height)

    Ngl.panel(wks,plots[n*2:n*2+2],[1,2],pres)

  Ngl.frame(wks)
  Ngl.destroy(wks)
  return(d)


if __name__ == "__main__":
  cscratch = os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton')
  ctl_dir =  os.path.join(cscratch, 'WRFv3.9_Sensitivity',
                          'WRFv3.9_Sensitivity_Ctl', 'WRFV3',
                          'run', 'summen_sensitivity_ctl')
  dry_dir = os.path.join(cscratch, 'WRFv3.9_Sensitivity',
                         'WRFv3.9_Sensitivity_DrySoil', 'WRFV3',
                         'run', 'summen_sensitivity_drysoil')
  vd = var_diff(os.path.join(ctl_dir, 'wrfsees_ccs3pb1_ls2_d02_2009-06*'),
                os.path.join(dry_dir, 'wrfsees_ccs3pb1_ls2_d02_2009-06*'),
                label_A = 'control',
                label_B = 'dry',
                varname='LH')
  read_data = True
  if read_data:
    vd.read_files()
    # TODO: this is a crude sychonization of the two time series.
    # Should generalize this.
    vd.data['control'] = vd.data['control'][12:273, ...]
    vd.time = vd.time[12:273]
  for this_t in range(2, 255):
    d = graphics(vd, t_idx=this_t, layer=None)
  Ngl.end()
