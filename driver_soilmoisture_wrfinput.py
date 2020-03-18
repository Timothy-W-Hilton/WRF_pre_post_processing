from adjust_WRF_soil_moisture import adjust_soil_moisture, gather_files
import os

wrf_run_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                           'twhilton', 'WRFv3.9_Sensitivity',
                           'FromScratch', 'WRFV3_HalfDry', 'run')
wrf_input_files = gather_files(wrf_run_dir, "wrfinput_d*")

for this_file in wrf_input_files:
    adjust_soil_moisture(this_file, 0.7,
                         soil_moist_vars=("SMOIS", ))
