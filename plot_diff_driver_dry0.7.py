import os
from plot_diff import var_diff, graphics

if __name__ == "__main__":
    rootdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                           'twhilton', 'WRFv3.9_Sensitivity', 'FromScratch')
    dry0p05_dir = os.path.join(rootdir, 'WRFV3_Dry', 'run')
    dry0p70_dir = os.path.join(rootdir, 'WRFV3_0.7Dry', 'run')
    ctl_dir = os.path.join(rootdir, 'WRFV3_Ctl', 'run')

    vd_0p05_0p70 = var_diff(os.path.join(dry0p05_dir, 'summen_sensitivity_dry',
                                         '*d01_2009-06-0[12]*'),
                            os.path.join(dry0p70_dir,
                                         'summen_sensitivity_0.7dry',
                                         '*d01_2009-06-0[12]*'),
                            label_A='dry 0.05',
                            label_B='dry 0.70',
                            varname='QCLOUD')

    vd_d0p7_ctl = var_diff(os.path.join(dry0p70_dir,
                                        'summen_sensitivity_0.7dry',
                                        '*d02_2009-06-0[12]*'),
                           os.path.join(ctl_dir,
                                        'summen_sensitivity_ctl',
                                        '*d02_2009-06-0[12]*'),
                           label_A='dry 0.70',
                           label_B='control',
                           varname='QCLOUD')
    read_data = True
    if read_data:
        vd_0p05_0p70.read_files()
        vd_d0p7_ctl.read_files()
        # vd_0p05_0p70.mask_land_or_water(mask_water=True)
    for this_t in range(0, 95):  # 255
        fig = graphics(vd_0p05_0p70, t_idx=this_t, layer=0, domain=2,
                       pfx='dry_0p05_0p70_diff')
        fig2 = graphics(vd_d0p7_ctl, t_idx=this_t, layer=0, domain=2,
                       pfx='dry_0p70_ctl_diff')
