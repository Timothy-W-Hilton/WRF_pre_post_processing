"""try out some visualizations of fog change/urbanization

uses universal tranverse mercator (UTM) zones to project the
coordinates.  This allows me to use the vectorized distance
calculation in geopandas geoseries (which, in turn, uses shapely's
distance function).  Shapely calculates cartesian distances, not
geodesic distances.  So the cost of using the vectorized distance
function is whatever loss of accuracy using UTM zones and cartesian
distances introduces.  For the purpose of sorting fog reduction by
distance to coast these distances look plenty close enough to me.
"""

import numpy as np
import pandas as pd
import seaborn as sns
from cartopy.feature import GSHHSFeature

from matplotlib.colors import cnames, Normalize
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io import shapereader
import shapely
import geopandas as gpd


def dist_to_coast_demo(na_w_coast, wrf_pixels):
    """plot a map showing WRF pixels colored by distance to coast

    Demonstrates that distance to coast calculation is working well
    enough.
    """
    fig = plt.figure()
    ax = plt.subplot(111)
    na_w_coast.plot(ax=ax)
    cm = ax.scatter(wrf_pixels['lon'],
                    wrf_pixels['lat'],
                    c=wrf_pixels['d_coast_km'],
                    s=4)
    cbar = plt.colorbar(cm)
    ax.set_xlabel('longitude ($^\circ$W)')
    ax.set_ylabel('latitude ($^\circ$N)')
    cbar.ax.set_title('km to coast')
    return(fig, ax)


def dfog_durban_scatterplot(df, show_ci=None, show_fits=False):
    """draw a scatterplot showing d(fog fraction) vs d(urban fraction)

    Each point is a WRF pixel.

    Uses the color axis to show distance to coast in km.

    returns:
    matplotlib figure and axes objects containing the plot
    """

    # draw the plot
    # fig = plt.figure()
    # ax = plt.subplot(111)
    n = Normalize().autoscale(A=df['d_coast_km'])
    bins = np.array([0, 5, 10, 15, 20, 25, 50, 100, 300, 500, 700])
    bins = np.array([0, 2, 5, 10, 15, 20, 25, 800])
    df['d_coast_binned'] = pd.cut(
        df['d_coast_km'],
        bins=bins,
    )

    tijuana = 32.522499
    santa_barbara = 34.4215357
    santa_cruz = 36.974474
    santa_rosa = 38.444660
    lat_bins = np.array([tijuana, santa_barbara, santa_cruz, santa_rosa,
                         41, 43, 45, 47, 49, 51])
    n_bins = lat_bins.size - 1
    df['lat_bins'] = pd.cut(df['lat'], lat_bins)

    sp = sns.lmplot(y="d_fog",
                    x="d_urban_frac",
                    hue='lat_bins',
                    palette=sns.color_palette("cubehelix",   # Blues_d",
                                              n_colors=n_bins),
                    data=df,
                    legend_out=True,
                    ci=show_ci,
                    fit_reg=show_fits)
    sp._legend.set_title('latitude $^\circ$N')
    fig = plt.gcf()
    ax = plt.gca()
    ax.set_xlim((-0.02, 1.02))
    ax.set_ylabel('$\Delta$fraction of hours with fog')
    ax.set_xlabel('$\Delta$urban fraction')
    return(fig, ax)

if __name__ == "__main__":

    df = pd.read_csv('./fog_change_data_frame.csv')

    # Read the Natural Earth shapefile dataset
    #----------------------------------
    kw = dict(resolution='10m', category='physical', name='coastline')
    coast_shp = shapereader.natural_earth(**kw)
    shp = shapereader.Reader(coast_shp)

    # define the coordinate reference system to latitude/longitude
    crs_latlon = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    crs_utmz10 = "+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    crs_utmz11 = "+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    coasts = gpd.read_file(coast_shp)
    coasts.crs = crs_latlon
    wrf_pixels = gpd.GeoDataFrame(
        df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=crs_latlon)
    wrf_pixels_utm = wrf_pixels.to_crs(crs_utmz10)
    # window coasts to W coast of N America
    na_w_coast = coasts.cx[-140:-110, 20:60]


    df =pd.DataFrame(
        {'City': ['San Diego', 'San Francisco', 'Truckee'],
         'Country': ['USA', 'USA', 'USA'],
         'Latitude': [32.715736, 37.733795, 39.327962],
         'Longitude': [-117.161087, -122.446747, -120.183253]})
    gdf = gpd.GeoDataFrame(df,
                           geometry=gpd.points_from_xy(df.Longitude,
                                                       df.Latitude),
                           crs=crs_latlon)
    ca_cities = gdf.to_crs(crs_utmz10)
    na_w_coast_utm = na_w_coast.to_crs(crs_utmz10)
    # window to California coast for testing
    ca_coast = na_w_coast_utm.iloc[15:17, :]
    # d1 = ca_cities.distance(ca_coast)
    # ax = ca_coast.plot()
    # ca_cities.plot(color='red', marker='x', ax=ax)

    # geopandas distance method operates on a *GeoSeries* object.  Not
    # a DataFrame, not a point, only a GeoSeries.  The geometry field
    # of a geopandas dataframe is a GeoSeries.
    for idx, this_city in ca_cities.iterrows():
        d = ca_coast.geometry.distance(this_city.geometry).values.min()
        d_km = d / 1000.0
        if d.size != 1:
            raise(ValueError('distance() returned multiple values'))
        print('{} to coast: {:0.0f} km'.format(this_city['City'], d_km))


    for idx, this_pixel in wrf_pixels.iterrows():
        d = ca_coast.geometry.distance(this_pixel.geometry).values
        d_km = 1 / 1000.0

    d = [ca_coast.geometry.distance(this_pixel.geometry).values.min() for
         idx, this_pixel in wrf_pixels_utm.iterrows()]
    d_km = np.array(d) / 1000.0
    wrf_pixels['d_coast_km'] = d_km

    # fig1, ax1 = dist_to_coast_demo(na_w_coast, wrf_pixels)
    wrf_pixels_coastal = wrf_pixels[wrf_pixels['d_coast_km'] < 5]
    fig2, ax2 = dfog_durban_scatterplot(wrf_pixels, show_fits=True)
    ax2.set_title('all WRF pixels with some urban use', pad=1)
    fig2.savefig('/Users/tim/work/Plots/Summen/NoUrban/deurbanize_fraction_vs_fog_change_allurban_latbins_lmfits.pdf')
