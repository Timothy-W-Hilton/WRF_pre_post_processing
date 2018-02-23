import os
from plot_diff import var_diff, graphics

if __name__ == "__main__":
    rootdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                           'twhilton', 'WRFv3.9_Sensitivity', 'FromScratch')
    ctl_dir = os.path.join(rootdir, 'WRFV3_Ctl', 'run')
    ctl2_dir = os.path.join(rootdir, 'WRFV3_Ctl2', 'run')

    vd = var_diff(os.path.join(ctl_dir, 'summen_sensitivity_ctl',
                               '*d01_2009-06*'),
                  os.path.join(ctl2_dir, 'summen_sensitivity_ctl2',
                               '*d01_2009-06*'),
                  label_A='ctl',
                  label_B='ctl2',
                  varname='QCLOUD')
    read_data = True
    if read_data:
        vd.read_files()
        # vd.mask_land_or_water(mask_water=True)
    for this_t in range(0, 120):  # 255
        fig = graphics(vd, t_idx=this_t, layer=0, domain=1,
                       pfx='ctl2')
