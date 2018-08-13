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

    def interpolate(self, data, method='NN', time_ax=0):
        """interpolate data

        Interpolate data using one of the supported methods.
        Supported methods include nearest neighbor ("NN") and inverse
        distance weighted ("IWD").

        ARGS
        method (str): {"NN"}|"IDW" Interploation method to use.  Must
           be one of "NN" (default) or "IDW".
        data (array-like): data to interpolate
        time_ax (int): axis (zero-based) containing the time dimension
           (default is 0)
        """
        self.inds = None
        if method is "NN":
            self._get_NN_idx()
        elif method is "IDW":
            raise NotImplementedError("IDW not yet implemented")
        else:
            raise ValueError("method argument must be one of ('NN', 'IWD')")
        idx2D = np.unravel_index(self.inds, dims=data[0, ...].shape)
        return(data[:, idx2D[0], idx2D[1]])
