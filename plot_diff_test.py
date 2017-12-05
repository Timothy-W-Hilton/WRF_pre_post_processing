import os
from glob import glob

import plot_diff as pld

if __name__ == "__main__":
    datadir = os.path.join('/', 'global', 'cscratch1', 'sd',
                           'twhilton', 'WRFv3.9_Sensitivity',
                           'FromScratch', 'WRFV3_urban', 'run',
                           'summen_sensitivity_urban')

    files = glob(os.path.join(datadir, "*d02*"))
    smois = pld.wrf_var(files, label='test', varname='SMOIS', is_atm=False)
    uvmet = pld.wrf_var(files, label='test2', varname='uvmet', is_atm=False)
    dzs = pld.wrf_var(files, label='test2', varname='DZS', is_atm=False)
    smois.read_files()
    uvmet.read_files()
    # dzs.read_files()    # fails due to bug - https://github.com/NCAR/wrf-python/issues/37
    smois.read_files(mask_land=True)
    smois.read_files(mask_water=True)
