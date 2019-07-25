"""change for Yatir forest WRF grid cells to Yatir land surface parameters

I implemented the Yatir Forest parameterization from Yosef et al
(2018) table S2 as land type 21 in the MODIFIED_IGBP_MODIS_NOAH land
use scheme for NOAH-MP.  This code finds WRF cells containing the
Yatir Forest eddy covariance tower () and for those cells:
    (1) changes the dominant land use to 21
    (2) changes the landuse fraction for type 21 to 1.0 and the
        fraction for all other land use types to 0.0


"""

import os
import warnings
from adjust_WRF_inputs import use_yatir_parameterization, gather_files

if __name__ == "__main__":
    wrf_run_dir = os.path.join(
        '/', 'global', 'cscratch1', 'sd',
        'twhilton', 'WRFv4.1_Experiments',
        'WRFv4.1_yatir_2015Aug_yatirparams_Z050_NCEPFNL_NOAHMP_ndom3',
        'WRFV4', 'run')

    wrfinput_files = sorted(gather_files(wrf_run_dir, "wrfinput*"))
    wrfbdy_files = sorted(gather_files(wrf_run_dir, "wrfbdy*"))
    wrflow_files = sorted(gather_files(wrf_run_dir, "wrflow*"))
    met_em_files_d02 = sorted(gather_files(wrf_run_dir, "met_em.d02*"))
    met_em_files_d03 = sorted(gather_files(wrf_run_dir, "met_em.d03*"))

    # for this_file in met_em_files_d01:
    #     print("starting met_em")
    #     use_yatir_parameterization(this_file)
    # for this_file in met_em_files_d02:
        # use_yatir_parameterization(this_file)
    for this_file in met_em_files_d03:
        use_yatir_parameterization(this_file)
    print("starting wrfinput")
    # warnings.warn("make sure to adjust wrfinput!", UserWarning)
    for this_file in wrfinput_files:
        use_yatir_parameterization(this_file)
    for this_file in wrfbdy_files:
        use_yatir_parameterization(this_file)
    for this_file in wrflow_files:
        use_yatir_parameterization(this_file)
