import sys
import numpy as np
import os
from scipy.stats import norm

from plot_diff import var_diff, MyFig
from map_tools_twh.map_tools_twh import CoastalSEES_WRF_prj
from map_tools_twh.map_tools_twh import CoastalSEES_WRF_Mapper
import matplotlib.cm
import matplotlib.patches as mpatches

DOMAIN = 2
read_vardiff_from_netcdf = True
adj_autocorr = True
significance_default = 0.50


def main(significance=None):

    if significance is None:
        significance = significance_default

    cscratchdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton')
    rootdir = os.path.join(cscratchdir, 'WRFv3.9_Sensitivity')
    ctl_dir = os.path.join(rootdir, 'WRFV3_Ctl',
                           'run', 'summen_sensitivity_ctl')
    urb_dir = os.path.join(rootdir, 'WRFV3_redwoodsurban',
                           'run', 'summen_redwoodsurban')
    # urb_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
    #                        'twhilton', 'WRFv3.9_Sensitivity',
    #                        'FromScratch', 'WRFV3_urban', 'run',
    #                        'summen_sensitivity_urban')
    out_dir = os.path.join(cscratchdir, 'plots_temporary')

    if read_vardiff_from_netcdf:
        vd = var_diff(ncfile=os.path.join('/', 'global', 'cscratch1', 'sd',
                                          'twhilton',
                                          'fog_pct_d02_RWurban.nc'))
    else:
        # wildcard_pat = "*d{:02d}02_2009-06-\{0[0-9],1[01]\}".format(DOMAIN)
        wildcard_pat = "*d{:02d}_2009-06-*".format(DOMAIN)
        vd = var_diff(os.path.join(ctl_dir, wildcard_pat),
                      os.path.join(urb_dir, wildcard_pat),
                      label_A='ctl',
                      label_B='urbanRW',
                      # varname='HFX')
                      varname='fogpresent')
                      # varname='fogbase')
                      # varname='fogpct')
        vd.read_files()
        vd.to_netcdf(os.path.join('/', 'global', 'cscratch1', 'sd',
                                  'twhilton',
                                  'fog_pct_d02_RWurban.nc'))

    print('calculating probabilities')
    z = vd.diff_means_test(adj_autocorr=adj_autocorr)
    vectorized_cdf = np.vectorize(lambda x: norm.cdf(x, 0.0, 1.0))
    p = vectorized_cdf(z)
    p = np.ma.masked_less(p, significance)
    # draw map
    print('drawing map')
    fig = MyFig(figsize=(8, 8))
    ax = fig.add_subplot(111,
                         projection=CoastalSEES_WRF_prj())
    this_cmap = matplotlib.cm.get_cmap('Dark2_r')
    this_map = CoastalSEES_WRF_Mapper(ax=ax, domain=DOMAIN)
    this_map.pcolormesh(vd.lon, vd.lat, p, vmin=0, vmax=1.0, cmap=this_cmap)
    ax.set_title(('urbanized redwood experiment\n'
                  'fog differences significant '
                  'at {:0.0f}%'.format(significance * 100.0)))
    significant_patch = mpatches.Patch(color=this_cmap.colors[-1],
                                       label='signficant')
    ax.legend(handles=(significant_patch, ))
    fname = ('fogpct_RWurban_d{domain:02d}_means_diff_'
             'map_2tail_adj{adj}_p{sig:0.2f}.png'.format(
                 domain=DOMAIN,
                 adj=adj_autocorr,
                 sig=significance))
    fig.savefig(fname=fname)
    print('saved {}'.format(fname))
    return(p)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        p = main(float(sys.argv[1]))
    else:
        p = main()
