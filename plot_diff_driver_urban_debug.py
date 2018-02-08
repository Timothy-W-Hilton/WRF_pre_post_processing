import os
from plot_diff import var_diff, VarDiffPlotter

DOMAIN = 2


def sum_layers(vd):
    """sum var_diff instance data for first 9 layers

    this is a quick and dirty way to look at cloud water content below
    400 m (the definition used by Johnstone and Dawson).  layer 10 is
    the first layer above 400 m in the Coastal SEES domain.
    """
    print('summing layers 0 to 9')
    for k, v in vd.data.items():
        vd.data[k] = v[:, 0:9, ...].sum(axis=1, keepdims=True)
        vd.data[k] = v.mean(axis=0, keepdims=True)
    vd.longname = vd.longname + ' below 400 m, time avg'
    return(vd)


if __name__ == "__main__":
    rootdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                           'twhilton', 'WRFv3.9_Sensitivity')
    ctl_dir = os.path.join(rootdir, 'WRFv3.9_Sensitivity_Ctl',
                           'WRFV3', 'run')
    ctl_dir_new = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton', 'WRFv3.9_Sensitivity',
                               'FromScratch', 'WRFV3_Ctl',
                               'run')
    urb_dir = os.path.join(rootdir, 'FromScratch', 'WRFV3_urban',
                           'run')

    vd = var_diff(os.path.join(ctl_dir_new, 'summen_sensitivity_ctl',
                               '*d{:02d}_2009-06-02*.nc'.format(DOMAIN)),
                  os.path.join(urb_dir, 'summen_sensitivity_urban',
                               '*d{:02d}_2009-06-02*.nc'.format(DOMAIN)),
                  label_A='ctl',
                  label_B='urban',
                  # varname='LCL')
                  varname='QCLOUD')
                  # varname='uvmet')
                  # varname='wa')
                  # varname='LU_INDEX')
    read_data = True
    if read_data:
        vd.read_files()
        # vd.mask_land_or_water(mask_water=False)
    for this_t in range(0, 47):  #
        plotter = VarDiffPlotter(vd, t_idx=this_t, layer=0,
                                 domain=DOMAIN,
                                 pfx='QCLOUD_debug04')
        if DOMAIN == 1:
            cb_orientation = 'horizontal'
        else:
            cb_orientation = 'vertical'
        fig = plotter.plot(cb_orientation=cb_orientation)

                           # vmin=0,
                           # vmax=500)
