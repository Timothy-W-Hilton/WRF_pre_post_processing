import numpy as np
import datetime
import os
from scipy.stats import norm

from plot_diff import var_diff, MyFig
from map_tools_twh.map_tools_twh import CoastalSEES_WRF_prj
from map_tools_twh.map_tools_twh import CoastalSEES_WRF_Mapper

DOMAIN = 2
read_data = True
adj_autocorr = True

if __name__ == "__main__":
    t0 = datetime.datetime.now()

    cscratchdir = os.path.join('/', 'global', 'cscratch1', 'sd',
                               'twhilton')
    rootdir = os.path.join(cscratchdir, 'WRFv3.9_Sensitivity')
    ctl_dir = os.path.join(rootdir, 'FromScratch', 'WRFV3_Ctl',
                           'run', 'summen_sensitivity_ctl')
    urb_dir = os.path.join(rootdir, 'FromScratch', 'WRFV3_redwoodsurban',
                           'run', 'summen_redwoodsurban')
    out_dir = os.path.join(cscratchdir, 'plots_temporary')

    if read_data:
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

    z = vd.diff_means_test(adj_autocorr=adj_autocorr)
    vectorized_cdf = np.vectorize(lambda x: norm.cdf(x, 0.0, 1.0))
    p = 1.0 - vectorized_cdf(z)

    # draw map
    fig = MyFig(figsize=(8, 8))
    ax = fig.add_subplot(111,
                         projection=CoastalSEES_WRF_prj())
    this_map = CoastalSEES_WRF_Mapper(ax=ax, domain=DOMAIN)
    this_map.pcolormesh(vd.lon, vd.lat, p, vmin=0, vmax=1.0)
    ax.set_title(r'Fog pct p(ctl > urban redwoods)')
    this_map.colorbar()
    fig.savefig(fname=('fogpct_d{domain:02d}_means_diff_'
                       'map_2tail_adj{adj}.png'.format(
                           domain=DOMAIN,
                           adj=adj_autocorr)))
