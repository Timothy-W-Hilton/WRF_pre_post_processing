# TODO:
# 1) make wrf_run_dir a command line argument
# 2) investigate why the previous job failed.  If timeout, use current
#    approach (adjust start time, resubmit).  If seg fault, lower
#    timestep and resubmit.
# 3) check last restart file to make sure it is complete.  If not,
# delete it and rerun get_last_restart_file()

import shutil
import wrf
import pandas as pd
import numpy as np
import netCDF4
import glob
import os
import f90nml
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

    def check_is_valid(self, fname_rst1, fname_rst0):
        """check whether a WRF restart file is "complete"

        Occasionally a WRF run times out in the middle of writing a
        restart file, leaving a valid-ish looking partial restart file
        on disk.  These incomplete restart files do not contain enough
        information to restart WRF, and cause WRF to crash.  It is
        therefore useful to identify these files when determininig at
        which timestamp WRF should restart.

        Tests whether the size of the restart file (in bytes) is
        within 2% of the previous restart file.  Turns out that
        checking that the length of the Times dimension was greater
        than 0 did not work; no sooner had I implemented this than WRF
        spit out a partial restart file with len(Times) of 1.

        During the Yatir runs, WRF is hitting a sementation fault and
        spitting out restart files that do not end on an even half
        hour interval.  Also check that that the timestamp on the
        restart file ends in 00 or 30.

        ARGS:
        fname_rst1 (string): full path of WRF restart file to be checked
        fname_rst0 (string): full path of WRF restart file preceding
            fname_rst1.

        RETURNS
        True if the restart file is complete, False otherwise

        """
        sz1 = os.path.getsize(fname_rst1)
        sz0 = os.path.getsize(fname_rst0)
        diff_ratio = abs((sz1 - sz0) / sz1)
        # consider fname_rst1 complete if size is within 2% of fname_rst0
        restart_file_is_large_enough = (diff_ratio < 0.02)
        if not restart_file_is_large_enough:
            print("{} is not large enough.".format(fname_rst1))
        # make sure fname_rst1 contains the correct 30-minute interval
        restart_file_is_even_half_hour = ((fname_rst1[-2:] == "00") or
                                          (fname_rst1[-2:] == "30"))
        if not restart_file_is_even_half_hour:
            print("{} is not an even 30-minutes.".format(fname_rst1))

        return(restart_file_is_large_enough and
               restart_file_is_even_half_hour)

    def get_last_restart_file(self, ndom, wildcard_str="wrfrst_d*[0-9][0-9]"):
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
            this_domain_str = 'd{:02}'.format(this_domain)
            restart_files[this_domain_str] = (
                sorted(glob.glob(os.path.join(self.run_dir,
                                              wildcard_str))))
            n_last_rst = min(map(len, restart_files.values()))
            last_restart_file = restart_files[this_domain_str][n_last_rst - 1]
            if self.check_is_valid(
                    last_restart_file,
                    restart_files[this_domain_str][n_last_rst - 2]):
                return(self.get_tstamp(last_restart_file))
            else:
                print(('{} is not valid for restarting WRF.'
                       '  Deleting and trying again'.format(
                           last_restart_file)))
                os.remove(last_restart_file)
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

    def get_new_start_time(self, new_start_time, verbose=True):
        """return a dict in containing the new start time

        The dict is formatted for insertion into a WRF namelist file
        using the f90nml module.

    """
        new_start_time = pd.to_datetime(new_start_time)
        namelist_update = {'time_control':
                           {'start_month': [new_start_time.month] * ndom,
                            'start_day': [new_start_time.day] * ndom,
                            'start_hour': [new_start_time.hour] * ndom,
                            'start_minute': [new_start_time.minute] * ndom,
                            'start_second': [new_start_time.second] * ndom}}
        if verbose:
            print("setting start time to {}".format(new_start_time))
        return(namelist_update)

    def update_namelist(self, namelist_update):
        """update values in a WRF namelist, save old namelist file as a backup

        The existing namelist file is copied to a new file with the
        current date appended to the file name (in format
        YYYY-MM-DD_hhmm)

        ARGS:
        namelist_update (dict): the new values to be placed in the
           namelist file.  Should be formatted as {category: {item:
           value}}.

        example:
        update_namelist({'domains': {'time_step": 100}})

        """
        fname_bak = self.fname + pd.Timestamp.now().strftime(".%Y-%m-%d_%H%M")
        shutil.copy(self.fname, fname_bak)
        f90nml.patch(fname_bak, namelist_update, self.fname)

    def get_ndomains(self):
        """return number of domains
        """
        return(self.nml['domains']['max_dom'])

    def get_current_timestep(self):
        """return time step from the namelist files
        """
        return(self.nml['domains']['time_step'])

    def get_new_timestep(self,
                         wrf_run_dir,
                         parent_job_id=None,
                         default_time_step=120):
        """update the WRF time_step value

        If the previous run ended in a segmentation fault, the time
        step needs to be shortened to try to satisfy the CFL condition
        (see,
        e.g. https://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy).
        As currently implemented, the time step is reduced to 2/3 its
        current value for another try.

        The previous run stdout/stderr file is assumed to be in a file
        named slurm-{parent_job_id}.out within the WRF run directory.
        Segmentation faults are identified by the presence of the
        string "Exited with exit code 174" in this file.

        If the previous run did not end in a segmentation fault or no
        previous run is provided returns a caller-specified default
        value (by default 120).  This is useful because shorter
        time_step values cause WRF to use more CPU time, so we want
        the maximum value that satisfies the CFL condition.

        ARGS:
        wrf_run_dir (string): full path to the WRF run directory
        parent_job_id (int): the SLURM job number for the previous run
        default_time_step (int): time_step value to use if previous
           run was successful (default is 120)

        RETURNS:
        dict containing the new timestep value in a format to be
        inserted into the WRF namelist file using f90nml.

        """
        time_step = self.get_current_timestep()
        parent_ended_in_segfault = False
        if (parent_job_id is not None) and (parent_job_id != 'None'):
            # if a parent file was provided, check its stdout/stderr
            # file to see if it ended in a segmentation fault.  If so,
            # reduce the time step for the next run.
            parent_outfile = os.path.join(wrf_run_dir,
                                          "slurm-{}.out".format(parent_job_id))
            for line in open(parent_outfile, 'r'):
                parent_ended_in_segfault = "Exited with exit code 174" in line
                if parent_ended_in_segfault:
                    break
        if parent_ended_in_segfault:
            # if segmentation fault happened, reduce the time step
            time_step = np.int(np.floor(time_step * 0.667))
        else:
            time_step = default_time_step

        print("setting time step to {}".format(time_step))
        return({'domains': {'time_step': time_step}})


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
    parser.add_argument('--default_wrf_timestep',
                        dest='default_time_step',
                        action='store',
                        default=120,
                        type=int,
                        help=(('default value for WRF timestep (seconds)'
                               '(default is 120)')))
    parser.add_argument('--parent_job_id',
                        dest='parent_job_id',
                        action='store',
                        default='None',
                        help=(('SLURM job id of the job that submitted the current task.  Used to '
                               'build the SLURM stdout/stderr filename of that job to determine why '
                               'the job "failed".  If failure was caused by timeout, then update the '
                               'start time and resubmit.  If failure was caused by segmentation fault, '
                               'then update the timestep and start time and resubmit.  ')))
    args = parser.parse_args()

    print("DEBUG PYTHON: parent_job_id: {}".format(args.parent_job_id))

    nml = WRF_namelist_file_tools(os.path.join(args.wrf_run_dir, args.fname))
    ndom = nml.get_ndomains()
    end_time = max(nml.get_end_time())

    rst = WRF_restart_files(args.wrf_run_dir)
    new_start_time = rst.get_last_restart_file(ndom)
    namelist_updates = {}
    if (new_start_time < end_time):
        print("updating {file} to start at {time}".format(
            file=os.path.join(args.wrf_run_dir, args.fname),
            time=new_start_time))
        namelist_updates.update(nml.get_new_start_time(new_start_time))
    namelist_updates.update(nml.get_new_timestep(
        args.wrf_run_dir,
        args.parent_job_id,
        default_time_step=args.default_time_step))
    nml.update_namelist(namelist_updates)
