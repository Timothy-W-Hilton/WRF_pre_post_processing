"""modified example from wrf module documentation

modified from https://wrf-python.readthedocs.io/en/latest/plot.html#plotting-a-two-dimensional-field

projection parameters of the WRF domain:

print(get_cartopy(slp).proj4_params)
{'a': 6370000.0, 'b': 6370000.0, 'nadgrids': '@null', 'proj': 'lcc', 'lon_0': -127.5, 'lat_0': 42.0, 'x_0': 0.0, 'y_0': 0.0, 'lat_1': 30.0, 'lat_2': 60.0}
"""

import os
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, ll_to_xy)


SanFrancisco = (37.7707405, -122.450006)  # lat deg N, lon deg E
SantaCruz = (36.974474, -122.028986)


def prepare_figure(smooth_slp):
    fig = plt.figure()
    ax = plt.axes(projection=get_cartopy(smooth_slp))
    states = NaturalEarthFeature(category="cultural", scale="10m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('10m', linewidth=0.8)
    ax.set_xlim(cartopy_xlim(smooth_slp))
    ax.set_ylim(cartopy_ylim(smooth_slp))
    # box containing San Francisco Bay, Monterey Bay
    ax.set_extent((-123.0, -121.0, 36.2, 38.0))
    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")
    return(fig, ax)


def plot_markers(locations, ax):
    """plot X for Santa Cruz, San Francisco

    locations: iterable of (lat, lon) tuples
    """
    for here in locations:
        ax.scatter(*reversed(here),
                   transform=crs.PlateCarree(),
                   marker='x',
                   c='blue')


# Open the NetCDF file
ncfile = Dataset(os.path.join(
    '/', 'Users', 'tim', 'work', 'Data',
    'SummenWRF',
    'metem_2005_ctl_NCEPDOE_d02_2005-06-01_00:00:00.nc'))

# Get the sea level pressure
slp = getvar(ncfile, "slp")

# Smooth the sea level pressure since it tends to be noisy near the
# mountains
smooth_slp = smooth2d(slp, 3, cenweight=4)

# Get the latitude and longitude points
lats, lons = latlon_coords(slp)

# Get the cartopy mapping object
cart_proj = get_cartopy(slp)

# Make the contour outlines and filled contours for the smoothed sea level
# pressure.
# plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp), 10, colors="blue",
#             transform=crs.PlateCarree())
# plt.contourf(to_np(lons), to_np(lats), to_np(smooth_slp), 10,
#              transform=crs.PlateCarree(),
#              cmap=get_cmap("Blues"))

# Add a color bar
# plt.colorbar(ax=ax, shrink=.98)

# plot markers for San Francisco, Santa Cruz

dummy = np.zeros(slp.shape)
for here in [SanFrancisco, SantaCruz]:
    xu, yu = to_np(ll_to_xy(ncfile, *here, stagger='U'))
    xv, yv = to_np(ll_to_xy(ncfile, *here, stagger='V'))
    x, y = to_np(ll_to_xy(ncfile, *here, stagger=None))
    print((xu, yu, xv, yv, x, y))
    dummy[yv, xv] = 500
# dummy = ma.masked_less(dummy, 500)

mylons = (lons[:, :-1] - (np.diff(lons, axis=1) / 2))[:-1, :]
mylats = (lats[:-1, :] - (np.diff(lats, axis=0) / 2))[:, :-1]

fig, ax = prepare_figure(smooth_slp)
ax.imshow(dummy,
          alpha=0.5,
          extent=(*cartopy_xlim(slp),
                  *cartopy_ylim(slp)),
          transform=get_cartopy(slp),
          origin='lower')
plot_markers((SanFrancisco, SantaCruz), ax)
plt.title("imshow")

fig, ax = prepare_figure(smooth_slp)
plt.pcolormesh(to_np(mylons),
               to_np(mylats),
               dummy[:-2, :-2],
               transform=crs.PlateCarree(),
               alpha=0.5,
               edgecolor='black')
plot_markers((SanFrancisco, SantaCruz), ax)
plt.title("pcolormesh, 1/2 lon/lat diffs substracted from mass grid coords")

plt.show()

# # prove that 'x' marker is plotted with center at x, y
# plt.figure()
# ax = plt.axes()
# ax.scatter(range(10), range(10), marker='x', s=1000)
# d = 0.02
# ax.set_xlim((3 - d, 3 + d))
# ax.set_ylim((3 - d, 3 + d))
