import os
from plot_diff import var_diff, graphics

if __name__ == "__main__":
    rootdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                           'twhilton', 'WRFv3.9_Sensitivity', 'FromScratch')
    ctl_dir = os.path.join(rootdir, 'WRFV3_Ctl', 'run')
    dry_dir = os.path.join(rootdir, 'WRFV3_Dry', 'run')
    # vd = var_diff(os.path.join(ctl_dir, 'met_em.d01.2009-06-02_00:00:00.nc'),
    #               os.path.join(dry_dir, 'met_em.d01.2009-06-02_00:00:00.nc'),
    #               label_A='ctl met_em',
    #               label_B='dry met_em',
    #               varname='SM000010')

    vd = var_diff(os.path.join(ctl_dir, 'summen_sensitivity_ctl',
                               'wrfsees_ccs3pb1_ls2_d02_2009-06-02_00:00:00'),
                  os.path.join(dry_dir, 'summen_sensitivity_dry',
                               'wrfsees_ccs3pb1_ls2_d02_2009-06-02_00:00:00'),
                  label_A='ctl out',
                  label_B='dry out',
                  varname='SMOIS')
    read_data = True
    if read_data:
        vd.read_files()
        # vd.mask_land_or_water(mask_water=True)
    for this_t in range(0, 6):  # 255
        fig = graphics(vd, t_idx=this_t, layer=0, domain=2)
