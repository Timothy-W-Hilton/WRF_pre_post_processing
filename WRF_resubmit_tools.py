# TODO:
# 1) make wrf_run_dir a command line argument
# 2) investigate why the previous job failed.  If timeout, use current approach (adjust start time, resubmit).  If seg fault, lower timestep and resubmit.
# 3) check last restart file to make sure it is complete.  If not,
# delete it and rerun get_last_restart_file()

import shutil
import wrf
import pandas as pd
import netCDF4
import glob
import os
import f90nml
import subprocess

wrf_run_dir = os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton',
                           'WRFv4.0_Sensitivity', 'WRFCLMv4.0_NCEPDOEp2_deurbanized',
                           'WRFV4', 'run')
fname_wrf_slurm = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton', 'WRFv4.0_Sensitivity',
                               'WRFCLMv4.0_wrfpy_test', 'WRFV4',
                               'run', 'WRF_Ctl_wrf.slurm')
fname = "namelist.input"


class WRF_restart_files(object):
    """tools to work with WRF restart files

    ARGS:
    fname (str): full path to WRF restart file
    """
    def __init__(self, run_dir):
        self.run_dir = run_dir

    def get_tstamp(self, fname):
    """read the timestamp from a WRF restart file

    assumes that the restart file contains a single timestamp in the
    netCDF variable Times.  If this assumption is violated results are
    undefined.

    ARGS:
    fname (str): full path to WRF restart file

    RETURNS:
    numpy.Datetime64 object containing the restart file timestamp
    """
    nc = netCDF4.Dataset(fname, 'r')
        self.tstamp = pd.to_datetime(wrf.extract_times(nc, 0))
    nc.close()
        return(self.tstamp)

    def get_last_restart_file(self, ndom, wildcard_str="wrfrst_d*"):
    """search a directory for most recent WRF restart file

    return the last time stamp that has a restart file for all
    domains.

    ARGS:
    ndom (int > 0): number of domains to search
    wildcard_str (str): wild card string matching restart file names
       (default "wrfrst_d*")

    RETURNS:
    numpy.Datetime64 object containing timestamp of last restart file
    """
    restart_files = {}
        # create a dict containing a list of WRF restart files present
        # for every WRF domain
    for this_domain in range(1, ndom + 1):
        restart_files['d{:02}'.format(this_domain)] = (
                sorted(glob.glob(os.path.join(self.run_dir,
                                              wildcard_str))))
    n_last_rst = min(map(len, restart_files.values()))
    last_restart_file = restart_files['d01'][n_last_rst - 1]
        return(self.get_tstamp(last_restart_file))


class WRF_namelist_file_tools(object):
    """tools for reading, updating WRF namelist files
    """
    def __init__(self, fname):
        """construct a WRF_namelist_file_tools object

        ARGS:
        fname: full path to the WRF namelist file
        """
        self.fname = fname
        self.nml = f90nml.read(fname)

    def get_end_time(self):
        """parse the WRF run end time to a pandas.Timestamp object

        RETURNS:
        pandas.Timestamp object containing the WRF run end time
        """
        end_time = map(pd.Timestamp,
                       self.nml['time_control']['end_year'],
                       self.nml['time_control']['end_month'],
                       self.nml['time_control']['end_day'],
                       self.nml['time_control']['end_hour'],
                       self.nml['time_control']['end_minute'],
                       self.nml['time_control']['end_second'])
        return(end_time)

    def update_namelist_start_time(self, new_start_time):
    """update the start time in a WRF namelist to a new value

    ARGS:
    namelist_file (string): full path to WRF namelist file to be updated
    new_start_time (numpy.Datetime64): the new start time
    """
    new_start_time = pd.to_datetime(new_start_time)
    namelist_update = {'time_control':
                       {'start_month': [new_start_time.month] * ndom,
                        'start_day': [new_start_time.day] * ndom,
                        'start_hour': [new_start_time.hour] * ndom,
                        'start_minute': [new_start_time.minute] * ndom,
                        'start_second': [new_start_time.second] * ndom}}
        fname_bak = self.fname + pd.Timestamp.now().strftime(".%Y-%m-%d_%H%M")
        shutil.copy(self.fname, fname_bak)
        f90nml.patch(fname_bak, namelist_update, self.fname)

    def get_ndomains(self):
        """return number of domains
        """
        return(self.nml['domains']['max_dom'])


if __name__ == "__main__":
    nml = WRF_namelist_file_tools(os.path.join(wrf_run_dir, fname))
    ndom = nml.get_ndomains()
    end_time = max(nml.get_end_time())

    rst = WRF_restart_files(wrf_run_dir)
    new_start_time = rst.get_last_restart_file(ndom)
    if (new_start_time < end_time):

    print("updating {file} to start at {time}".format(
            file=os.path.join(wrf_run_dir, fname), time=new_start_time))
    nml.update_namelist_start_time(new_start_time)
        # print('running `sbatch {slurm_script}` now'.format(
        #     slurm_script=fname_wrf_slurm))
        # subprocess.run(["sbatch", fname_wrf_slurm], stdout=subprocess.PIPE)
