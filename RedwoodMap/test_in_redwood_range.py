import os.path
import fiona
from shapely.geometry import shape,mapping, Point, Polygon, MultiPolygon

fname = os.path.join('/', 'Users', 'tim', 'work', 'Data',
                     'Redwoods',
                     'Redwood_SequoiaSempervirens_extentNorthAmerica',
                     'data', 'commondata', 'data0', 'sequsemp.shp')
rr_polys = fiona.open(fname)  # Redwood range polygons
rr_shapes = [shape(this_poly['geometry']) for this_poly in rr_polys]

bbrwsp = Point(37.1737, -122.2464)  # Big Basin Redwoods State Park



for i, this_shape in enumerate(rr_shapes):
    in_poly = bbrwsp.within(this_shape)
    print("{:02d} {}".format(i, in_poly))
