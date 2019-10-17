import pandas as pd
from io import StringIO
from tabulate import tabulate

vegparm_str = """
MODIFIED_IGBP_MODIS_NOAH   ALBD,   SLMO,   SFEM,   SFZ0, THERIN,   SCFX,   SFHC, LANDCOVER
1,     12.,   .30,   .95,   50.,    4.,  3.33, 29.2e5, 'Evergreen Needleleaf Forest'
2,     12.,   .50,   .95,   50.,    5.,  1.67, 29.2e5, 'Evergreen Broadleaf Forest'
3,     14.,   .30,   .94,   50.,    4.,  2.86, 25.0e5, 'Deciduous Needleleaf Forest'
4,     16.,   .30,   .93,   50.,    4.,  2.63, 25.0e5, 'Deciduous Broadleaf Forest'
5,     13.,   .30,   .97,   50.,    4.,  2.11, 41.8e5, 'Mixed Forests'
6,     22.,   .10,   .93,    5.,    3.,  1.56, 20.8e5, 'Closed Shrublands'
7,     20.,   .15,   .95,    6.,    3.,  2.14, 20.8e5, 'Open Shrublands'
8,     22.,   .10,   .93,    5.,    3.,  1.56, 20.8e5, 'Woody Savannas'
9,     20.,   .15,   .92,   15.,    3.,  2.00, 25.0e5, 'Savannas'
10,    19.,   .15,   .96,   12.,    3.,  2.37, 20.8e5, 'Grasslands'
11,    14.,   .42,   .95,   30.,   5.5,  1.32, 35.5e5, 'Permanent wetlands'
12,    17.,   .30,  .985,   15.,    4.,  2.71, 25.0e5, 'Croplands'
13,    15.,   .10,   .88,   80.,    3.,  1.67, 18.9e5, 'Urban and Built-Up'
14,    18.,   .25,   .98,   14.,    4.,  2.56, 25.0e5, 'cropland/natural vegetation mosaic'
15,    55.,   .95,   .95,   0.1,    5.,    0., 9.0e25, 'Snow and Ice'
16,    25.,   .02,   .90,    1.,    2.,  0.81, 12.0e5, 'Barren or Sparsely Vegetated'
17,     8.,   1.0,   .98,  0.01,    6.,    0., 9.0e25, 'Water'
18,    15.,   .50,   .93,   30.,    5.,  2.67, 9.0e25, 'Wooded Tundra'
19,    15.,   .50,   .92,   15.,    5.,  2.67, 9.0e25, 'Mixed Tundra'
20,    12.,   .02,   .97,   50.,     3.,  1.67, 18.9e5, 'Yatir'
21,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
22,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
23,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
24,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
25,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
26,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
27,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
28,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
29,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
30,    15.,   .02,   .88,   80.,    3.,  1.67, 18.9e5, 'Unassigned'
31,    10.,   .10,   .97,   80.,    3.,  1.67, 18.9e5, 'Low Intensity Residential '
32,    10.,   .10,   .97,   80.,    3.,  1.67, 18.9e5, 'High Intensity Residential'
33,    10.,   .10,   .97,   80.,    3.,  1.67, 18.9e5, 'Industrial or Commercial'
"""

landparm_str = """MODIFIED_IGBP_MODIS_NOAH
SHDFAC, NROOT,   RS,      RGL,      HS,      SNUP,  MAXALB,   LAIMIN,  LAIMAX,   EMISSMIN, EMISSMAX, ALBEDOMIN, ALBEDOMAX,   Z0MIN,   Z0MAX,   ZTOPV,    ZBOTV, LANDCOVER'
1,       .70,   4,    125.,    30.,   47.35,   0.08,    52.,    5.00,   6.40,   .950,    .950,     .12,      .12,      .50,     .50,      17.0,     8.5,    'Evergreen Needleleaf Forest'
2,      .95,   4,    150.,    30.,   41.69,   0.08,    35.,    3.08,   6.48,   .950,    .950,     .12,      .12,      .50,     .50,      35.0,     1.0,    'Evergreen Broadleaf Forest'
3,      .70,   4,    150.,    30.,   47.35,   0.08,    54.,    1.00,   5.16,   .930,    .940,     .14,      .15,      .50,     .50,      14.0,     7.0,    'Deciduous Needleleaf Forest'
4,      .80,   4,    100.,    30.,   54.53,   0.08,    58.,    1.85,   3.31,   .930,    .930,     .16,      .17,      .50,     .50,      20.0,    11.5,    'Deciduous Broadleaf Forest'
5,      .80,   4,    125.,    30.,   51.93,   0.08,    53.,    2.80,   5.50,   .930,    .970,     .17,      .25,      .20,     .50,      18.0,    10.0,    'Mixed Forests'
6,      .70,   3,    300.,   100.,   42.00,   0.03,    60.,    0.50,   3.66,   .930,    .930,     .25,      .30,      .01,     .05,      0.50,    0.10,    'Closed Shrublands'
7,      .70,   3,    170.,   100.,   39.18,  0.035,    65.,    0.60,   2.60,   .930,    .950,     .22,      .30,      .01,     .06,      0.50,    0.10,    'Open Shrublands'
8,      .70,   3,    300.,   100.,   42.00,   0.03,    60.,    0.50,   3.66,   .930,    .930,     .25,      .30,      .01,     .05,      0.50,    0.10,    'Woody Savannas'
9,      .50,   3,     70.,    65.,   54.53,   0.04,    50.,    0.50,   3.66,   .920,    .920,     .20,      .20,      .15,     .15,      0.50,    0.10,    'Savannas'
10,     .80,   3,     40.,   100.,   36.35,   0.04,    70.,    0.52,   2.90,   .920,    .960,     .19,      .23,      .10,     .12,      0.50,    0.01,    'Grasslands'
11,      .60,   2,     70.,    65.,   55.97,   0.015,     59.,    1.75,   5.72,   .950,    .950,     .14,      .14,      .30,     .30,      0.00,    0.00,    'Permanent wetlands'
12,     .80,   3,     40.,   100.,   36.25,   0.04,    66.,    1.56,   5.68,   .920,    .985,     .17,      .23,      .05,     .15,      0.50,    0.01,    'Croplands'
13,     .10,   1,    200.,   999.,   999.0,   0.04,    46.,    1.00,   1.00,   .880,    .880,     .15,      .15,      .50,     .50,      0.00,    0.00,    'Urban and Built-Up'
14,      .80,   3,     40.,   100.,   36.25,   0.04,    68.,    2.29,   4.29,   .920,    .980,     .18,      .23,      .05,     .14,      0.50,    0.01,    'cropland/natural vegetation mosaic'
15,     .00,   1,    999.,   999.,   999.0,   0.02,    82.,    0.01,   0.01,   .950,    .950,     .55,      .70,    0.001,   0.001,      0.00,    0.00,    'Snow and Ice'
16,     .01,   1,    999.,   999.,   999.0,   0.02,    75.,    0.10,   0.75,   .900,    .900,     .38,      .38,      .01,     .01,      0.02,    0.01,    'Barren or Sparsely Vegetated'
17,     .00,   0,    100.,    30.,   51.75,   0.01,    70.,    0.01,   0.01,   .980,    .980,     .08,      .08,   0.0001,  0.0001,      0.00,    0.00,    'Water'
18,     .60,   3,    150.,   100.,   42.00,  0.025,    55.,    0.41,   3.35,   .930,    .930,     .15,      .20,      .30,     .30,      10.0,    0.10,    'Wooded Tundra'
19,     .60,   3,    150.,   100.,   42.00,  0.025,    60.,    0.41,   3.35,   .920,    .920,     .15,      .20,      .15,     .15,      5.00,    0.10,    'Mixed Tundra'
20,     .56,   4,    286.,    30.,   47.35,   0.08,    52.,    3.00,   3.00,   .970,    .970,     .12,      .12,      2.0,     2.0,      10.0,     1.0,    'Yatir'
"""

def format_vegparm():
    vegparm = pd.read_csv(header=1, sep=',', filepath_or_buffer=StringIO(vegparm_str))
    vegparm_formatted = tabulate(vegparm, headers='keys', tablefmt='pipe')
    return(vegparm_formatted)

def format_landparm():
    landparm = pd.read_csv(header=1, sep=',', filepath_or_buffer=StringIO(landparm_str))
    landparm_formatted = tabulate(landparm, headers='keys', tablefmt='pipe')
    return(landparm_formatted)
