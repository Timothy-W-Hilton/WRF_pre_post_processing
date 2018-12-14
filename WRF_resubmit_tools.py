import wrf
import pandas as pd
import netCDF4
import glob
import os
import f90nml

wrf_run_dir = os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton',
                           'WRFv4.0_Sensitivity', 'WRFCLMv4.0_wrfpy_test',
                           'WRFV4', 'run')
# wrf_run_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
#                            'twhilton', 'WRFv4.0_Sensitivity',
#                            'WRFCLMv4.0_NCEPDOEp2_deurbanized',
#                            'WRFV4', 'run', '')
fname = "namelist.input"


def get_restart_file_tstamp(fname):
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
    tstamp = wrf.extract_times(nc, 0)
    nc.close()
    return(tstamp)


def get_last_restart_file(wrf_run_dir, ndom, wildcard_str="wrfrst_d*"):
    """search a directory for most recent WRF restart file

    return the last time stamp that has a restart file for all
    domains.

    ARGS:
    wrf_run_dir (str): full path to the WRF run directory to be
       searched for restart files
    ndom (int > 0): number of domains to search
    wildcard_str (str): wild card string matching restart file names
       (default "wrfrst_d*")

    RETURNS:
    numpy.Datetime64 object containing timestamp of last restart file
    """
    restart_files = {}
    for this_domain in range(1, ndom + 1):
        restart_files['d{:02}'.format(this_domain)] = (
            sorted(glob.glob(os.path.join(wrf_run_dir, wildcard_str))))
    n_last_rst = min(map(len, restart_files.values()))
    last_restart_file = restart_files['d01'][n_last_rst - 1]
    return(get_restart_file_tstamp(last_restart_file))


def update_namelist_start_time(namelist_file, new_start_time):
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
    f90nml.patch(namelist_file, namelist_update, namelist_file + ".new")


if __name__ == "__main__":
    nml = f90nml.read(os.path.join(wrf_run_dir, fname))
    ndom = nml['domains']['max_dom']

    # f90nml.patch(os.path.join(wrf_run_dir, fname),
    #              new_start_time,
    #              os.path.join(wrf_run_dir, "namelist_test.input"))
    new_start_time = get_last_restart_file(wrf_run_dir, ndom)
    print("updating {file} to start at {time}".format(
        file=fname, time=new_start_time))
    update_namelist_start_time(os.path.join(wrf_run_dir, fname),
                               str(new_start_time))
