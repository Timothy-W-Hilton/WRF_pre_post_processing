import os
from glob import glob
import netCDF4
from wrf import getvar, ALL_TIMES

if __name__ == "__main__":

    datadir = os.path.join('/', 'global', 'cscratch1',
                           'sd', 'twhilton', 'WRFv3.9_Sensitivity',
                           'WRFv3.9_Sensitivity_Ctl', 'WRFV3', 'run',
                           'summen_sensitivity_ctl')
    files = glob(os.path.join(datadir, "wrfsees_ccs3pb1_ls2_d02_*"))
    nclist = [netCDF4.Dataset(f, mode="r") for f in files]
    smois = getvar(nclist, 'SMOIS', timeidx=ALL_TIMES)
    smois = getvar(nclist, 'SMOIS', timeidx=ALL_TIMES)
    uvmet = getvar(nclist, 'uvmet', timeidx=ALL_TIMES)
    uvmet10 = getvar(nclist, 'uvmet10', timeidx=ALL_TIMES)
    ua = getvar(nclist, 'ua', timeidx=ALL_TIMES)
    va = getvar(nclist, 'va', timeidx=ALL_TIMES)
    lat = getvar(nclist, 'lat', timeidx=ALL_TIMES)
    lon = getvar(nclist, 'lon', timeidx=ALL_TIMES)
    zs = getvar(nclist, 'ZS', timeidx=ALL_TIMES, meta=False)[0, ...]
    dzs = getvar(nclist, 'DZS', timeidx=ALL_TIMES, meta=False)[0, ...]
    island = getvar(nclist, 'XLAND', timeidx=ALL_TIMES, meta=False)[0, ...]
    for this_nc in nclist:
        this_nc.close()
