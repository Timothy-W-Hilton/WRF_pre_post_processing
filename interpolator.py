import numpy as np
from scipy.spatial import cKDTree


class cKDTreeInterpolator(object):

    def __init__(self, lon_in, lat_in, lon_out, lat_out):
        """
        """
        self.lon_in = lon_in
        self.lat_in = lat_in
        self.lon_out = lon_out
        self.lat_out = lat_out
        self._build_tree()

    def _lon_lat_to_cartesian(self, lon, lat, R=1):
        """
        Convert spherical coordinates to three-dimensional Cartesian
        coordinates.

        calculates three dimensional cartesian coordinates (x, y, z) for
        specified longitude, latitude coordinates on a sphere with radius
        R.  Written and posted at earthpy.org by Oleksandr Huziy.
        http://earthpy.org/interpolation_between_grids_with_ckdtree.html
        accessed 19 Mar 2014 by Timothy W. Hilton.

        PARAMETERS
        ==========
        lon; np.ndarray: longitude values
        lat; np.ndarray: latitude values
        R: scalar; radius of the sphere.

        RETURNS
        =======
        three element tuple containing X, Y, and Z, one element per
        lon,lat pair.
        """
        lon_r = np.radians(lon)
        lat_r = np.radians(lat)

        x = R * np.cos(lat_r) * np.cos(lon_r)
        y = R * np.cos(lat_r) * np.sin(lon_r)
        z = R * np.sin(lat_r)
        return x, y, z

    def _build_tree(self):
        """
        find nearest neighbors for a set of lon, lat points from a second
        set of lon, lat points.
        """

        # convert spherical lon, lat coordinates to cartesian coords. Note
        # that these x,y,z are 3-dimensional cartesian coordinates of
        # positions on a sphere, and are thus different from the x,y,z
        # *indices* of the STEM grid.
        self.xi, self.yi, self.zi = self._lon_lat_to_cartesian(self.lon_in,
                                                               self.lat_in)
        # use a K-dimensional tree to find the nearest neighbor to (x,y,z)
        # from the points within (xs, ys, zs).  A KD tree is a data
        # structure that allows efficient queries of K-dimensional space (K
        # here is 3).
        self.tree = cKDTree(np.dstack((self.xi.flatten(),
                                       self.yi.flatten(),
                                       self.zi.flatten())).squeeze())

    def _get_NN_idx(self):
        """
        """
        self.xo, self.yo, self.zo = self._lon_lat_to_cartesian(self.lon_out,
                                                               self.lat_out)
        d, inds = self.tree.query(
            np.dstack((self.xo, self.yo, self.zo)).squeeze(), k=1)
        # self.inds = np.unravel_index(inds, self.lon_out.shape)
        self.inds = inds

    def interpolate(self, data, method='NN'):
        """interpolate data

        Interpolate data using one of the supported methods.
        Supported methods include nearest neighbor ("NN") and inverse
        distance weighted ("IWD").

        ARGS
        method (str): {"NN"}|"IDW" Interploation method to use.  Must
           be one of "NN" (default) or "IDW".
        data (array-like): data to interpolate
        """
        self.inds = None
        if method is "NN":
            self._get_NN_idx()
        elif method is "IDW":
            raise NotImplementedError("IDW not yet implemented")
        else:
            raise ValueError("method argument must be one of ('NN', 'IWD')")

        #brute force array filling
        new_data = np.full(self.inds.shape, np.nan)
        # f = open('coords.txt', 'w')
        for i in range(self.inds.shape[0]):
            for j in range(self.inds.shape[1]):
                new_data[i, j] = data.flatten()[self.inds[i, j]]
                thisidx = self.inds[i, j]
        #         print(thisidx,
        #               self.lon_in.flatten()[thisidx],
        #               self.lat_in.flatten()[thisidx],
        #               self.lon_out[i, j],
        #               self.lat_out[i, j],
        #               data.flatten()[thisidx],
        #               new_data[i, j],
        #               file=f)
        # f.close()
        # self.debug_plot()
        return(new_data)
        # return(data.flatten()[self.inds])

    def debug_plot(self):
        """
        """
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        fig, ax = plt.subplots(
            ncols=2,
            nrows=2,
            subplot_kw={'projection': ccrs.PlateCarree()})
        ax[0, 0].scatter(self.lon_in, self.lat_in, marker='.')
        ax[0, 0].set_title('grid in')
        ax[0, 0].add_feature(cfeature.LAND)
        ax[0, 0].add_feature(cfeature.OCEAN)

        ax[0, 1].scatter(self.lon_out, self.lat_out, marker='.')
        ax[0, 1].set_title('grid out')
        ax[0, 1].add_feature(cfeature.LAND)
        ax[0, 1].add_feature(cfeature.OCEAN)


def find_nearest_xy(lon_in, lat_in, lon_out, lat_out):
    """
    find nearest neighbors for a set of lon, lat points from a second
    set of lon, lat points.

    Given a set of arbitrary (lon, lat) positions, find the horizontal
    (x, y) STEM grid indices of the nearest STEM grid cell center to
    each position.
    PARAMETERS
    ----------
    lon, lat: ndarray; of arbitrary longitudes and latitudes.  Must
       contain the same number of elements.
    lon_stem, lat_stem: ndarrays; longitudes and latitudes of STEM
       grid cell centers. Must contain the same number of elements.

    RETURNS:
    an N-element tuple of X and Y indices, one index per observation.
        The indices are the closest point in (lon, lat) to each point
        in (lon_stem, lat_stem).  N is therefore equal to the number
        of dimensions in lon_stem and lat_stem.
    """
    # convert spherical lon, lat coordinates to cartesian coords. Note
    # that these x,y,z are 3-dimensional cartesian coordinates of
    # positions on a sphere, and are thus different from the x,y,z
    # *indices* of the STEM grid.
    lon_in, lat_in = np.meshgrid(lon_in, lat_in)
    lat_in = lat_in[::-1]
    xi, yi, zi = _lon_lat_to_cartesian(lon_in, lat_in)
    xo, yo, zo = _lon_lat_to_cartesian(lon_out, lat_out)

    # use a K-dimensional tree to find the nearest neighbor to (x,y,z)
    # from the points within (xs, ys, zs).  A KD tree is a data
    # structure that allows efficient queries of K-dimensional space (K
    # here is 3).
    tree = cKDTree(np.dstack((xi.flatten(),
                              yi.flatten(),
                              zi.flatten())).squeeze())

    d, inds = tree.query(
        np.dstack((xo, yo, zo)).squeeze(), k=1)

    return(np.unravel_index(inds, lon_in.shape))


def _lon_lat_to_cartesian(lon, lat, R=1):
    """
    Convert spherical coordinates to three-dimensional Cartesian
    coordinates.

    calculates three dimensional cartesian coordinates (x, y, z) for
    specified longitude, latitude coordinates on a sphere with radius
    R.  Written and posted at earthpy.org by Oleksandr Huziy.
    http://earthpy.org/interpolation_between_grids_with_ckdtree.html
    accessed 19 Mar 2014 by Timothy W. Hilton.

    PARAMETERS
    ==========
    lon; np.ndarray: longitude values
    lat; np.ndarray: latitude values
    R: scalar; radius of the sphere.

    RETURNS
    =======
    three element tuple containing X, Y, and Z, one element per
    lon,lat pair.
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x, y, z
