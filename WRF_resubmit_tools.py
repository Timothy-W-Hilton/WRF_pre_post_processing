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
import argparse

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

    def check_is_valid(self, fname_rst):
        """check whether a WRF restart file is "complete"

        Occasionally a WRF run times out in the middle of writing a
        restart file, leaving a valid-ish looking partial restart file
        on disk.  These incomplete restart files do not contain enough
        information to restart WRF, and cause WRF to crash.  It is
        therefore useful to identify these files when determininig at
        which timestamp WRF should restart.

        Currently tests whether the length of the "Times" netcdf
        dimension and variable is 1.  In my experience partial restart
        files have a Times dimension of 0.

        ARGS:
        fname_rst (string): full path of WRF restart file

        RETURNS
        True if the restart file is complete, False otherwise
        """
        nc = netCDF4.Dataset(fname_rst)
        ntimes = len(nc.variables['Times'])
        nc.close()
        return(ntimes == 1)

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
        if self.check_is_valid(last_restart_file):
        return(self.get_tstamp(last_restart_file))
        else:
            shutil.delete(last_restart_file)
            return(self.get_last_restart_file(ndom, wildcard_str))


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

    parser = argparse.ArgumentParser(description=(
        ("script to update a WRF namelist.input file to begin a new "
         "run where the previous run left off.")))
    parser.add_argument('--wrf_run_dir',
                        dest='wrf_run_dir',
                        action='store',
                        default='./',
                        help=('full path to the WRF run directory'))
    parser.add_argument('--namelist_input',
                        dest='fname',
                        action='store',
                        default='namelist.input',
                        help=(('name of the namelist file to be updated'
                               '(default is namelist.input)')))
    parser.add_argument('--parent_job_id',
                        dest='parent_job_id',
                        action='store',
                        default='',
                        help=(('SLURM job id of the job that submitted the current task.  Used to '
                               'build the SLURM stdout/stderr filename of that job to determine why '
                               'the job "failed".  If failure was caused by timeout, then update the '
                               'start time and resubmit.  If failure was caused by segmentation fault, '
                               'then update the timestep and start time and resubmit.  ')))
    args = parser.parse_args()

    nml = WRF_namelist_file_tools(os.path.join(args.wrf_run_dir, args.fname))
    ndom = nml.get_ndomains()
    end_time = max(nml.get_end_time())

    rst = WRF_restart_files(args.wrf_run_dir)
    new_start_time = rst.get_last_restart_file(ndom)
    if (new_start_time < end_time):

    print("updating {file} to start at {time}".format(
            file=os.path.join(args.wrf_run_dir, args.fname), time=new_start_time))
    nml.update_namelist_start_time(new_start_time)
        # print('running `sbatch {slurm_script}` now'.format(
        #     slurm_script=fname_wrf_slurm))
        # subprocess.run(["sbatch", fname_wrf_slurm], stdout=subprocess.PIPE)
