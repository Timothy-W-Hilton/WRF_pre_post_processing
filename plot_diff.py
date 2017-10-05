#
#  File:
#    newcolor4.py
#
#  Synopsis:
#    Illustrates new color capabilities in PyNGL 1.5.0
#
#  Categories:
#    contour plots
#
#  Author:
#    Mary Haley, based on NCL example
#
#  Date of initial publication:
#    November 2012
#
#  Description:
#    This example shows how to draw nine plots with three different colormaps
#
#  Effects illustrated:
#    o  Using the new "cnFillPalette" resource.
#    o  Panelling plots.
#    o  Drawing an Aitoff map
#    o  Generating dummy data
#
#  Output:
#     A single visualization with nine plots
#
import numpy as np
import os
import Ngl, Nio
import netCDF4

class var_diff(object):
  def __init__(self, f_tim, f_osu, varname):
    self.fnames = {'tim':f_tim, 'OSU':f_osu}
    self.varname = varname
    self.longname = None
    self.units = None
    self.lat = None
    self.lon = None
    self.data = {'tim':None, 'OSU':None}

  def read_files(self):
    """read variable from Tim output, OSU output
    """
    for k, v in self.data.iteritems():
      nf = netCDF4.Dataset(self.fnames[k])
      self.data[k] = nf[self.varname][...]
      # read latitude
      if self.lat is None:
        self.lat = nf['XLAT'][0, ...]
      else:
        if np.allclose(nf['XLAT'][0, ...],
                       self.lat,
                       equal_nan=True) is False:
          raise RuntimeError('OSU latitude differs from Tim latitude')
      # read longitude
      if self.lon is None:
        self.lon = nf['XLONG'][0, ...]
      else:
        if np.allclose(nf['XLONG'][0, ...],
                       self.lon,
                       equal_nan=True) is False:
          raise RuntimeError('OSU longitude differs from Tim longitude')
      # read units
      if self.units is None:
        self.units = nf[self.varname].units
      else:
        if nf[self.varname].units != self.units:
          raise RuntimeError('OSU units differs from Tim units')
      self.longname = nf[self.varname].description
      nf.close()

def graphics(vd):
  """draw contour maps of two different model runs' values, along with
  the difference

  ARGS:
  vd (var_diff object): values for one variable from two different model runs

  """
  #---Start the graphics section
  wks_type = "png"
  wks = Ngl.open_wks(wks_type,"{}_diff_maps".format(vd.varname))

  # One color map per row of plots
  colormaps = ["wgne15","StepSeq25","BlueDarkRed18"]

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
  res.pmTickMarkDisplayMode  = "Never"

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
  plots  = []

  all_data = np.stack( (vd.data['tim'], vd.data['OSU']) ).flatten()
  dmin = min(all_data)
  dmax = max(all_data)

  for k in vd.data.keys():
    print("plot {} data".format(k))

    #---This is a new resource added in PyNGL 1.5.0
    res.cnFillPalette          = colormaps[0]

    print("contour min: {}, max: {}".format(dmin, dmax))
    res.cnMinLevelValF         = dmin
    res.cnMaxLevelValF         = dmax
    res.cnMaxLevelCount        = 10

    plots.append(Ngl.contour_map(wks,vd.data[k][1, ...],res))

  # Resources for panelling
  pres                  = Ngl.Resources()
  pres.nglFrame         = False
  pres.nglPanelLabelBar = True

  # Calculate start Y position for first row of plots
  height = 0.45            # we know this will be height of small plots
  extra  = 1.0-(2*height)
  top    = 1.0-(extra/2.)

  # Draw a title before we draw plots
  title               = "Multiple panels on one page, 3 different colormaps"
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
  Ngl.end()

if __name__ == "__main__":
  vd = var_diff('tim_wrf3.8.1_sees_ccs3pb1_ls2_d01_2009-03-01_00_00_00',
                './OSU_wrfsees_ccs2_d01_2009-03-01_00_00_00',
                varname='LH')
  vd.read_files()
  graphics(vd)
