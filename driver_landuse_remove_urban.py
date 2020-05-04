"""reduce soil moisture in met_em* to 0.05 of off-the-shelf level

alters land use in wrfinput* files in
/global/cscratch1/sd/twhilton/WRFv3.9_Sensitivity/FromScratch/WRFV3_urban/run
"""

import os
import warnings
from adjust_WRF_inputs import remove_urban, gather_files

if __name__ == "__main__":
    # wrf_run_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
    #                            'twhilton', 'WRFv3.9_Sensitivity',
    #                            'FromScratch', 'WRFV3_nourbanNOAH', 'run')
    wrf_run_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton', 'WRFv4.0_Sensitivity',
                               'WRFCLMv4.0_NCEPDOEp2_deurbanized',
                               'WRFV4', 'run')
    # wrf_run_dir = "/global/cscratch1/sd/twhilton/nourban_test/"
    wrfinput_files = sorted(gather_files(wrf_run_dir, "wrfinput*"))
    met_em_files_d01 = sorted(gather_files(wrf_run_dir, "met_em.d01*"))
    met_em_files_d02 = sorted(gather_files(wrf_run_dir, "met_em.d02*"))

    # for this_file in met_em_files_d01:
    #     print("starting met_em")
    #     remove_urban(this_file)
    # for this_file in met_em_files_d02:
    #     print("starting met_em")
    #     remove_urban(this_file)
    print("starting wrfinput")
    # warnings.warn("make sure to adjust wrfinput!", UserWarning)
    remove_urban(wrfinput_files[0])
    remove_urban(wrfinput_files[1])
