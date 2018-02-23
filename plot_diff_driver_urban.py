import os
from plot_diff import var_diff, graphics

if __name__ == "__main__":
    rootdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                           'twhilton', 'WRFv3.9_Sensitivity', 'FromScratch')
    ctl_dir = os.path.join(rootdir, 'WRFV3_Ctl', 'run')
    urb_dir = os.path.join(rootdir, 'WRFV3_urban', 'run')
    # vd = var_diff(os.path.join(ctl_dir, 'met_em.d01.2009-06-02_00:00:00.nc'),
    #               os.path.join(urb_dir, 'met_em.d01.2009-06-02_00:00:00.nc'),
    #               label_A='ctl met_em',
    #               label_B='dry met_em',
    #               varname='SM000010')

    vd = var_diff(os.path.join(ctl_dir, 'summen_sensitivity_ctl',
                               '*d01_2009-06*'),
                  os.path.join(urb_dir, 'summen_sensitivity_urban',
                               '*d01_2009-06*'),
                  label_A='ctl out',
                  label_B='urban out',
                  varname='QCLOUD')
    read_data = True
    if read_data:
        vd.read_files()
        # vd.mask_land_or_water(mask_water=True)
    for this_t in range(20, 120):  # 255
        fig = graphics(vd, t_idx=this_t, layer=0, domain=1,
                       pfx='urban')
