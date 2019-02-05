"""modified example from wrf module documentation

modified from https://wrf-python.readthedocs.io/en/latest/plot.html#plotting-a-two-dimensional-field
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

# Create a figure
fig = plt.figure(figsize=(12,6))
# Set the GeoAxes to the projection used by WRF
ax = plt.axes(projection=cart_proj)

# Download and add the states and coastlines
states = NaturalEarthFeature(category="cultural", scale="10m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.coastlines('10m', linewidth=0.8)

# Make the contour outlines and filled contours for the smoothed sea level
# pressure.
# plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp), 10, colors="blue",
#             transform=crs.PlateCarree())
# plt.contourf(to_np(lons), to_np(lats), to_np(smooth_slp), 10,
#              transform=crs.PlateCarree(),
#              cmap=get_cmap("Blues"))

# Add a color bar
# plt.colorbar(ax=ax, shrink=.98)

# Set the map bounds
ax.set_xlim(cartopy_xlim(smooth_slp))
ax.set_ylim(cartopy_ylim(smooth_slp))
ax.set_extent((-123.0, -121.0, 36.2, 38.0))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title("Sea Level Pressure (hPa)")

# plot markers for San Francisco, Santa Cruz
SanFrancisco = (37.7707405, -122.450006)  # lat deg N, lon deg E
SantaCruz = (36.974474, -122.028986)
xy_coords = ll_to_xy(ncfile, *zip(*(SanFrancisco, SantaCruz)))
dummy = np.zeros(slp.shape)
for here in [SanFrancisco, SantaCruz]:
    ax.scatter(*reversed(here),
               transform=crs.PlateCarree(),
               marker='x',
               c='blue')
    dummy[tuple(to_np(ll_to_xy(ncfile, *here, stagger=None))[::-1])] = 500
dummy = ma.masked_less(dummy, 500)
plt.pcolormesh(to_np(lons),
               to_np(lats),
               dummy,
               transform=crs.PlateCarree(),
               edgecolor='black')
plt.show()
