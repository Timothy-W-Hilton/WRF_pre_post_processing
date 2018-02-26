"""process Redwood range digital data into shapely.geometry.Polygon objects
"""

import fiona
from shapely.geometry import shape
import geopandas as gp

def _redwood_range_to_shapes_list(fname_rw):
    """make list of shapely.geometry.shape objects from Redwoods range shapefile

    Helper function for redwood_range_to_gdf

    ARGS:
    fname_rw (str): full path to redwood location shapefile .shp file

    RETURNS:
    list of `shapely.geometry.shape` objects
    """
    rw_polys = fiona.open(fname_rw)  # Redwood range polygons
    rw_shapes = [shape(this_poly['geometry']) for this_poly in rw_polys]
    return(rw_shapes)

def redwood_range_to_gdf(fname_rw):
    """make geopandas.GeoDataFrame from Redwoods range shapefile

    ARGS:
    fname_rw (str): full path to redwood location shapefile .shp file

    RETURNS
    a `geopandas.GeoDataFrame` object containing the polygons
    comprising the redwood range
    """
    gdf = gp.GeoDataFrame(geometry=_redwood_range_to_shapes_list(fname_rw))
    return(gdf)
