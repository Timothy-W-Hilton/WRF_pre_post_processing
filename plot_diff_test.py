import os
from glob import glob

import plot_diff as pd

if __name__ == "__main__":
    datadir = os.path.join('/', 'global', 'cscratch1',
                           'sd', 'twhilton', 'WRFv3.9_Sensitivity',
                           'WRFv3.9_Sensitivity_Ctl', 'WRFV3', 'run',
                           'summen_sensitivity_ctl')
    files = glob(os.path.join(datadir, "wrfsees_ccs3pb1_ls2_d02_*"))
    v1 = pd.wrf_var(files, label='test', varname='SMOIS', is_atm=False)
    v1.read_files()
    # v1.read_files(mask_land=True)
    # v1.read_files(mask_water=True)
