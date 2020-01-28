"""delete most WRF restart files from a WRF run directory

I am writing out one ~150 MB WRF restart file per 30 simulation
minutes.  This enables immediate recovery from a segmentation fault
caused by a CFL condition violation (see
https://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition) in
the interest of using scarce CPU time efficiently.  It also uses a lot
of disk storage.

This module scans a specified WRF directory for WRF restart files.  It
keeps the ten most recent and one per day (at time 00:00:00), and
deletes the rest to preserve disk space.

WRF restart files are assumed to be named in the format
"wrfrst_d[0-9][0-9]*".

"""

import glob
import os
import argparse

restart_wildcard_string = "wrfrst_d{domain:02d}*"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=(
        ("Delete most WRF restart files from a WRF run directory.  The "
         "following files are kept: the ten most recent, and one restart "
         "file per day of the WRF run (with timestamp matching "
         "00:00).")))
    parser.add_argument('--wrf_run_dir',
                        dest='wrf_run_dir',
                        action='store',
                        default='./',
                        help=('full path to the WRF run directory'))
    parser.add_argument('--n_domains',
                        dest='ndom',
                        action='store',
                        default=2,
                        type=int,
                        help=(("number of nested domains in "
                              " the WRF run (default 2)")))
    args = parser.parse_args()

    for this_domain in range(1, args.ndom + 1):
        restart_files = sorted(glob.glob(
            os.path.join(args.wrf_run_dir,
                         restart_wildcard_string.format(domain=this_domain))))
        # keep the most recent ten files
        restart_files = restart_files[:-10]
        # keep files with timestamps of 00:00
        for this_file in restart_files:
            if "_00:00:00" not in this_file:
                print("deleting {}".format(os.path.basename(this_file)))
                os.remove(this_file)
