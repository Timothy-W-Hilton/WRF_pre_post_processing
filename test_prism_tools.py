import os.path
import prism_tools

if __name__ == "__main__":
    pts = prism_tools.PRISMTimeSeries(
        os.path.join('/', 'Users',
                     'tim', 'work', 'Data', 'PRISM',
                     'PRISM_tmean_stable_4kmD1_20090601_20090630_asc.zip'),
        varname='tmean', varunits='C')
    # pmp._find_data_files()
    # hdr_vars, data = pmp._parse_file(pmp.data_files[0])
    pts.from_zip_archive()
