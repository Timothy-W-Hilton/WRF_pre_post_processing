import os
from plot_diff import var_diff, graphics

if __name__ == "__main__":
    cscratch = os.path.join('/', 'global', 'cscratch1', 'sd', 'twhilton')
    ctl_dir = os.path.join(cscratch, 'WRFv3.9_Sensitivity',
                           'WRFv3.9_Sensitivity_Ctl_short', 'WRFV3',
                           'run')
    dry_dir = os.path.join(cscratch, 'WRFv3.9_Sensitivity',
                           'WRFv3.9_Sensitivity_DrySoil_newrst', 'WRFV3',
                           'run', 'summen_sensitivity_drysoil')
    vd = var_diff(os.path.join(ctl_dir, 'wrfrst_d01_2009-06-02*'),
                  os.path.join(dry_dir, 'wrfsees_ccs3pb1_ls2_d01_2009-06-02*'),
                  label_A='short_ctl rst',
                  label_B='short_dry out',
                  varname='QCLOUD')
    read_data = True
    if read_data:
        vd.read_files()
        # vd.mask_land_or_water(mask_water=True)
    for this_t in range(0, 1):  # 255
        fig = graphics(vd, t_idx=this_t, layer=0, domain=1)
