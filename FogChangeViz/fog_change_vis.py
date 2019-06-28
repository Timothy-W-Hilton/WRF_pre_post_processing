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

def dist_to_coast(df_row, shp=None):
    """calculate distance from a lat/lon point to a coastline
    """
    pt1 = shapely.geometry.Point(df_row['lon'], df_row['lat'])
    return(np.min([pt1.distance(this_g) for this_g in shp.geometries()]))

if __name__ == "__main__":

    df = pd.read_csv('./fog_change_data_frame.csv')
    df['round_lat'] = round(df['lat'])
    sns.scatterplot(y="d_fog", x="d_urban_frac", hue='round_lat', data=df)

    coastlines = GSHHSFeature()

    # Read the Natural Earth shapefile dataset
    #----------------------------------
    kw = dict(resolution='10m', category='physical', name='coastline')
    states_shp = shapereader.natural_earth(**kw)
    shp = shapereader.Reader(states_shp)
    pt1 = shapely.geometry.Point(df['lat'][0], df['lon'][0])
    # for this_g in shp.geometries():
    #     distance_between_pts = pt1.distance(this_g)
    #     print(distance_between_pts)
    d = np.min([pt1.distance(this_g) for this_g in shp.geometries()])
    # TODO: this needs projection information somehow...
    df['dist_to_coast'] = df.apply(dist_to_coast, axis=1, args={shp: shp})

    # draw the plot
    n = Normalize().autoscale(A=df['dist_to_coast'])
    df['dist_to_coast_cat'] = pd.cut(
        df['dist_to_coast'],
        bins=np.append(np.linspace(0.0, 1.0, 10), 10)
    )
    sns.scatterplot(y="d_fog",
                    x="d_urban_frac",
                    hue='dist_to_coast_cat',
                    hue_norm=n,
                    palette=sns.color_palette("Blues_d", n_colors=10),
                    data=df)
