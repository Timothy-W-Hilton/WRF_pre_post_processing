import os.path
import prism_tools
import pickle

parse_data = False

if __name__ == "__main__":
    prism_dir = os.path.join('/', 'Users',
                             'tim', 'work', 'Data', 'PRISM')
    fname_pickle = "pts.pickle"
    if parse_data:

        pts = prism_tools.PRISMTimeSeries(
            os.path.join(prism_dir,
                         'PRISM_tmean_stable_4kmD1_20090601_20090630_asc.zip'),
            varname='tmean', varunits='C')
        pts.from_zip_archive()
        pfile = open(fname_pickle, 'wb')
        pickle.dump(pts, pfile)
        pfile.close()
    else:
        pfile = open(fname_pickle, 'rb')
        pts = pickle.load(pfile)
        pfile.close()

    lon, lat = prism_tools.read_WRF_latlon(
        os.path.join(prism_dir, 'WRF_d02_latlon.nc'))
    # new_image = pts.interpolate(lon, lat)
