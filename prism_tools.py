"""parse PRISM data; resample to another grid

`PRISM <http://prism.oregonstate.edu>` is a elevation-sensitive
reconstruction of meteorological fields

"""
import re
import zipfile
import numpy as np

fp_tol = 1e-6   # floating point tolerance

class PRISMTimeSeries(object):
    """container for time series of `PRISM <http://prism.oregonstate.edu>`
    data with timestamps, coordinates

    """

    def __init__(self, fname):
        """intialize filename

        ARGS:
        fname: name of data file in zip archive
        """
        self.fname = fname
        self.lat = None
        self.lon = None
        self.tstamps = None
        self.data = None
        self.varname = None
        self.varunits = None


class PRISMMonthlyParser(object):
    """parse a zipped monthly `PRISM <http://prism.oregonstate.edu>`
    archive into a PRISMTimeSeries object

    handles reading of data files from `PRISM
    <http://prism.oregonstate.edu>` zipped archive.
    """
    def __init__(self, infile):
        """initialize filename

        ARGS:
        infile (str): full path to zip archive
        """
        self.infile = infile
        self.lat = None
        self.lon = None

    def _find_data_files(self):
        """locate data files ("*.asc") in a PRISM zip archive

        places list of data files in the archive into self.data_files
        """
        zfile = zipfile.ZipFile(self.infile, mode='r')
        files = zfile.namelist()
        zfile.close()
        re_datafile = re.compile('\.asc$')
        data_files = list(filter(re_datafile.search, files))
        self.data_files = data_files
        self.data_list = [None] * len(self.data_files)

    def _parse_file(self, fname, n_hdr_lines=6):
        """parse a daily PRISM data file from a zip archive
        """
        print("parsing {}".format(fname))
        zarchive = zipfile.ZipFile(self.infile, mode='r')
        zfile = zarchive.open(fname, mode='r')
        hdr_vars = {}
        hdr_lines = [next(zfile) for x in range(n_hdr_lines)]
        for line in hdr_lines:
            name, var = line.split()
            hdr_vars[name.strip()] = float(var)
        data = np.genfromtxt(zfile)
        zfile.close()
        zarchive.close()
        data[np.abs(data - hdr_vars[b'NODATA_value']) < 1e-6] = np.nan
        # make sure all the data were read
        assert(data.shape == (np.int(hdr_vars[b'nrows']),
                              np.int(hdr_vars[b'ncols'])))
        # calculate lat, lon coordinates
        lon_coords = hdr_vars[b'xllcorner'] + (np.arange(hdr_vars[b'nrows']) *
                                               hdr_vars[b'cellsize'])
        lat_coords = hdr_vars[b'yllcorner'] + (np.arange(hdr_vars[b'ncols']) *
                                               hdr_vars[b'cellsize'])
        if self.lat is None:
            self.lat = lat_coords
        else:
            assert(np.array_equal(self.lat, lat_coords))
        if self.lon is None:
            self.lon = lon_coords
        else:
            assert(np.array_equal(self.lon, lon_coords))
        return(data)

    def parse_all(self):
        """parse all ASCII data files from a PRISM monthly data zip archive
        """
        self._find_data_files()
        data_list = [self._parse_file(this_file) for
                     this_file in self.data_files]
        self.data = np.array(self.data_list)
