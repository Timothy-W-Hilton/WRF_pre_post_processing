import os
from plot_diff import var_diff, VarDiffPlotter

DOMAIN = 2


def is_foggy_obrien_2013(vd):
    """return true if any layer <= 400m has qc > 0.05 g / kg

    Implements the definition of fog from O'Brien et al 2013: fog is
    defintied is present in a horizontal gridcell if any layer at or
    below 400 m has liquid water content (qc) >= 0.05 g / kg.

    REFERENCES

    Brien, T. A., L. C. Sloan, P. Y. Chuang, I. C. Faloona, and
    J. A. Johnstone (2013), Multidecadal simulation of coastal fog
    with a regional climate model, Climate Dynamics, 40(11-12),
    2801-2812, doi:10.1007/s00382-012-1486-x.

    """
    vd.longname = 'fog_present'
    gram_per_kg = 1e-3
    vertical_axis = 1  # axes are (0=time, 1=vertical, 2=x, 3=y)
    for k, v in vd.data.items():
        layer_is_foggy = v[:, 0:9, ...] >= (0.05 * gram_per_kg)
        vd.data[k] = layer_is_foggy.any(axis=vertical_axis)
    return(vd)


def sum_layers(vd, time_avg=False):
    """sum var_diff instance data for first 9 layers

    this is a quick and dirty way to look at cloud water content below
    400 m (the definition used by Johnstone and Dawson).  layer 10 is
    the first layer above 400 m in the Coastal SEES domain.
    """
    print('summing layers 0 to 9')
    vd.longname = vd.longname + ' below 400 m'
    for k, v in vd.data.items():
        vd.data[k] = v[:, 0:9, ...].sum(axis=1, keepdims=True)
        if time_avg:
            vd.data[k] = v.mean(axis=0, keepdims=True)
    if time_avg:   # outside loop so string is only appended once
        vd.longname = vd.longname + ' time avg'
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
                               '*d{:02d}_2009-06*.nc'.format(DOMAIN)),
                  os.path.join(urb_dir, 'summen_sensitivity_urban',
                               '*d{:02d}_2009-06*.nc'.format(DOMAIN)),
                  label_A='ctl',
                  label_B='urban',
                  # varname='HFX')
                  # varname='LH')
                  # varname='LCL')
                  varname='QCLOUD')
                  # varname='uvmet')
                  # varname='wa')
                  # varname='LU_INDEX')
    read_data = True
    if read_data:
        vd.read_files()
        # vd.mask_land_or_water(mask_water=False)
    # vd = sum_layers(vd, time_avg=True)
    vd = is_foggy_obrien_2013(vd)
    pfx = 'urban'
    for this_series in ['all_tstamps', 'time_avg']:
        if this_series == 'all_tstamps':
            t_end = 1440
            pfx = '2018-02-09'
        else:
            t_end = 1
            pfx = pfx + '_timeavg'
            vd = sum_layers(vd, time_avg=True)
        for this_t in range(0, t_end):  #
            plotter = VarDiffPlotter(vd, t_idx=this_t, layer=0,
                                     domain=DOMAIN,
                                     pfx=pfx)
            if DOMAIN == 1:
                cb_orientation = 'horizontal'
            else:
                cb_orientation = 'vertical'
            fig = plotter.plot(cb_orientation=cb_orientation,
                               vmin=-0.0000001,
                               vmax=1.0000001)
