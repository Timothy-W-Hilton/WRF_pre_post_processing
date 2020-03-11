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
from adjust_WRF_inputs import use_yatir_parameterization
from adjust_WRF_inputs import replace_landusef_luindex
from adjust_WRF_inputs import gather_files
from adjust_WRF_inputs import adjust_soil_moisture

if __name__ == "__main__":
    # wrf_run_dir = os.path.join(
    #     '/', 'global', 'cscratch1', 'sd',
    #     'twhilton', 'WRFv4.1_Experiments',
    #     'WRFv4.1_yatir_2015Aug_yatirparams_NCEPFNL_NOAHMP_ndom3_wetsoil',
    #     'WRFV4', 'run')

    VWC_factor = 1.1  # adjustment factor for soil moisture
                      # (i.e. volumetric water content)

    # directory of files for testing
    wrf_run_dir = '/global/cscratch1/sd/twhilton/test_adjust_land/'

    wrfinput_files = sorted(gather_files(wrf_run_dir, "wrfinput*"))
    wrfbdy_files = sorted(gather_files(wrf_run_dir, "wrfbdy*"))
    wrflow_files = sorted(gather_files(wrf_run_dir, "wrflow*"))
    met_em_files_d01 = sorted(gather_files(wrf_run_dir, "met_em.d01*"))
    met_em_files_d02 = sorted(gather_files(wrf_run_dir, "met_em.d02*"))
    met_em_files_d03 = sorted(gather_files(wrf_run_dir, "met_em.d03*"))

    for this_file in met_em_files_d01:
        print("starting met_em")
        use_yatir_parameterization(this_file)
    for this_file in met_em_files_d02:
        use_yatir_parameterization(this_file)
        adjust_soil_moisture(this_file, VWC_factor, yatir_only=True)
    for this_file in met_em_files_d03:
        use_yatir_parameterization(this_file)
        adjust_soil_moisture(this_file, VWC_factor, yatir_only = True)
    print("starting wrfinput")
    # for this_wrfinput_file, this_metem_set in zip(wrfinput_files,
    #                                               (met_em_files_d01,
    #                                                met_em_files_d02,
    #                                                met_em_files_d03)):
    #     use_yatir_parameterization(this_wrfinput_file)
    #     replace_landusef_luindex(this_wrfinput_file,
    #                              this_metem_set[0])
    #     adjust_soil_moisture(this_wrfinput_file,
    #                          VWC_factor,
    #                          yatir_only=True,
    #                          soil_moist_vars=('SMOIS', ))
    # for this_file in wrfbdy_files:
    #     use_yatir_parameterization(this_file)
    # for this_file in wrflow_files:
    #     use_yatir_parameterization(this_file)
