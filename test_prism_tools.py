import os.path
import prism_tools

if __name__ == "__main__":
    prism_dir = os.path.join('/', 'Users',
                             'tim', 'work', 'Data', 'PRISM')
    pts = prism_tools.PRISMTimeSeries(
        os.path.join(prism_dir,
                     'PRISM_tmean_stable_4kmD1_20090601_20090630_asc.zip'),
        varname='tmean', varunits='C')
    pts.from_zip_archive()

    lat, lon = prism_tools.read_WRF_latlon(
        os.path.join(prism_dir, 'WRF_d02_latlon.nc'))
    pts.interpolate()
