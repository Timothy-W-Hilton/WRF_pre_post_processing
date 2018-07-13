import prism_tools

if __name__ == "__main__":
    pmp = prism_tools.PRISMMonthlyParser('/Users/tim/work/Data/PRISM/PRISM_tmean_stable_4kmD1_20090601_20090630_asc.zip')
    # pmp._find_data_files()
    # hdr_vars, data = pmp._parse_file(pmp.data_files[0])
    pmp.parse_all()
