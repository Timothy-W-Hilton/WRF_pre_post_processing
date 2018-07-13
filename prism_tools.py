"""parse PRISM data; resample to another grid

`PRISM <http://prism.oregonstate.edu>` is a elevation-sensitive
reconstruction of meteorological fields

"""
import re
import zipfile
import numpy as np
from datetime import datetime

fp_tol = 1e-6   # floating point tolerance

class PRISMTimeSeries(object):
    """container for time series of `PRISM <http://prism.oregonstate.edu>`
    data with timestamps, coordinates

    """

    def __init__(self, fname, varname, varunits):
        """intialize filename

        ARGS:
        fname: name of data file in zip archive
        """
        self.fname = fname
        self.lat = None
        self.lon = None
        self.tstamps = None
        self.data = None
        self.varname = varname
        self.varunits = varunits

    def from_zip_archive(self):
        # calculate lat, lon coordinates
        pmp = PRISMMonthlyParser(self.fname)
        pmp.parse_all()
        self.data = pmp.data
        self.tstamps = pmp.tstamp_list
        nrows = pmp.data.shape[1]
        ncols = pmp.data.shape[2]
        self.lon = pmp.xllcorner + (np.arange(nrows) * pmp.cellsize)
        self.lat = pmp.yllcorner + (np.arange(ncols) * pmp.cellsize)


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
        self.xllcorner = None
        self.yllcorner = None
        self.cellsize = None

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
        if self.xllcorner is None:
            self.xllcorner = hdr_vars[b'xllcorner']
        if self.yllcorner is None:
            self.yllcorner = hdr_vars[b'yllcorner']
        if self.cellsize is None:
            self.cellsize = hdr_vars[b'cellsize']
        # make sure all grid dimensions match
        assert(self.xllcorner == hdr_vars[b'xllcorner'])
        assert(self.yllcorner == hdr_vars[b'yllcorner'])
        assert(self.cellsize == hdr_vars[b'cellsize'])
        return(data)

    def parse_all(self):
        """parse all ASCII data files from a PRISM monthly data zip archive
        """
        self._find_data_files()
        data_list = [self._parse_file(this_file) for
                     this_file in self.data_files]
        self.data = np.array(data_list)
        self.tstamp_list = [self.get_time_stamps(this_file) for
                            this_file in self.data_files]

    def get_time_stamps(self, fname):
        """parse data from PRISM filename
        """
        re_tstamp = re.compile('[0-9]{8}')
        tstamp_str = re_tstamp.search(fname)
        import pdb; pdb.set_trace()
        if tstamp_str:
            tstamp = datetime.strptime(tstamp_str.group(0), '%Y%m%d')
            return(tstamp)
        else:
            warnings.warn('unable to parse time stamp from {}'.format(fname))
            return(None)
