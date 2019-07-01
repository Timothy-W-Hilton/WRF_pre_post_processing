"""try out some visualizations of fog change/urbanization
"""

import numpy as np
import pandas as pd
import seaborn as sns
from cartopy.feature import GSHHSFeature

from matplotlib.colors import cnames, Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io import shapereader
import shapely
import geopandas as gpd

def dist_to_coast(df_row, shp=None):
    """calculate distance from a lat/lon point to a coastline
    """
    pt1 = shapely.geometry.Point(df_row['lon'], df_row['lat'])
    return(np.min([pt1.distance(this_g) for this_g in shp.geometries()]))

if __name__ == "__main__":

    df = pd.read_csv('./fog_change_data_frame.csv')
    df['round_lat'] = round(df['lat'])

    # Read the Natural Earth shapefile dataset
    #----------------------------------
    kw = dict(resolution='10m', category='physical', name='coastline')
    coast_shp = shapereader.natural_earth(**kw)
    shp = shapereader.Reader(coast_shp)
    # TODO: this needs projection information somehow...
    # foo = df.copy()
    # df = df.head()
    # df['dist_to_coast'] = df.apply(dist_to_coast, axis=1, args={shp: shp})

    # # draw the plot
    # n = Normalize().autoscale(A=df['dist_to_coast'])
    # df['dist_to_coast_cat'] = pd.cut(
    #     df['dist_to_coast'],
    #     bins=np.append(np.linspace(0.0, 1.0, 10), 10)
    # )
    # sns.scatterplot(y="d_fog",
    #                 x="d_urban_frac",
    #                 hue='dist_to_coast_cat',
    #                 hue_norm=n,
    #                 palette=sns.color_palette("Blues_d", n_colors=10),
    #                 data=df)

    # 200-300 contains all of US west coast
    # search for lon < 50 and lat > 0 to (mostly) to (crudely) window N America

    # define the coordinate reference system to latitude/longitude
    crs_latlon = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    crs_utmz10 = "+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    coasts = gpd.read_file(coast_shp)
    coasts.crs = crs_latlon
    wrf_pixels = gpd.GeoDataFrame(
        df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=crs_latlon)
    # d1 = wrf_pixels.distance(coasts)
    # d2 = coasts.distance(wrf_pixels)
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
    ca_coast = na_w_coast_utm.iloc[16:17, :]
    d1 = ca_cities.distance(ca_coast)
    ax = ca_coast.plot()
    ca_cities.plot(color='red', marker='x', ax=ax)
    ca_cities.iloc[2].geometry.distance(ca_coast.iloc[0].geometry)
    ca_coast.iloc[0].distance(ca_cities.iloc[0])
    # geopandas distance method operates on a *GeoSeries* object.  Not
    # a DataFrame, not a point, only a GeoSeries.  The geometry field
    # of a geopandas dataframe is a GeoSeries.
    for idx, this_city in ca_cities.iterrows():
        d = ca_coast.geometry.distance(this_city.geometry).values
        d_km = d / 1000.0
        if d.size != 1:
            raise(ValueError('distance() returned multiple values'))
        print('{} to coast: {:0.0f} km'.format(this_city['City'], d_km[0]))


    foo = ca_coast.geometry.distance(ca_cities.iloc[0].geometry)
