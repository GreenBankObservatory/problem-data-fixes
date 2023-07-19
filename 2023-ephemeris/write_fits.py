from astropy.io import fits
import glob

def open_original_fits(fits_path):
    original_fits = fits.open(fits_path)
    return original_fits

def get_n_hdu(data_hdulist):
    return len(data_hdulist)

def check_hdr_key(k):
    keep_key = True

    keys_to_ignore = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 
                      'NAXIS2', 'PCOUNT', 'GCOUNT', 'TFIELDS']
    starts_to_ignore = ['TTYPE', 'TFORM', 'TUNIT']
    
    if k in keys_to_ignore:
        keep_key = False
    for s in starts_to_ignore:
        if k.startswith(s):
            keep_key = False
    
    return keep_key
    
def get_tbl_hdr(tbl):
    all_keys = {}
    for k in tbl.header.keys():
        if check_hdr_key(k):
            all_keys[k] = tbl.header[k]
    return all_keys

def write_tbl(data_hdulist, hdu_num, new_data, n_rows):
    old_tbl = data_hdulist[hdu_num]
    new_colnames = new_data.keys()
    new_tbl = fits.BinTableHDU.from_columns(old_tbl.columns, nrows=n_rows)
    new_hdr = get_tbl_hdr(old_tbl)
    for new_col in new_colnames:
        new_tbl.data[new_col] = new_data[new_col]
    for new_key in new_hdr.keys():
        new_tbl.header.set(new_key, new_hdr[new_key])
    data_hdulist[hdu_num] = new_tbl
    return data_hdulist

def write_hdr(data_hdulist, hdu_num, new_hdr):
    for k in new_hdr.keys():
        data_hdulist[hdu_num].header.set(k, new_hdr[k])
    return data_hdulist

def write_file(fits_file, save_path, overwrite_bool=False):
    fits_file.writeto(save_path, overwrite=overwrite_bool)

def write_new_LO1A(old_path, new_path, new_hdr0, new_data3, n_rows3):
    # OPEN THE ORIGINAL FILE
    lo1a_fits = open_original_fits(old_path)
    # (0) PRIMARY HDU'S HEADER
    lo1a_fits = write_hdr(lo1a_fits, 0, new_hdr0)
    # (3) LO1TBL HDU'S DATA
    lo1a_fits = write_tbl(lo1a_fits, 3, new_data3, n_rows3)
    # SAVE THE NEW FILE
    write_file(lo1a_fits, new_path, overwrite_bool=True)
    # CLOSE THE ORIGINAL FILE
    lo1a_fits.close()

def write_new_VEGAS(old_path, new_path, new_hdr0, new_data6, n_rows6):
    # new_hdr0, new_data1, n_rows1, new_hdr4, new_data4, n_rows4, 
    # OPEN THE ORIGINAL FILE
    vegas_fits = open_original_fits(old_path)
    # (0) PRIMARY HDU'S HEADER
    vegas_fits = write_hdr(vegas_fits, 0, new_hdr0)
    # (1) SPURS HDU'S DATA
    #vegas_fits = write_tbl(vegas_fits, 1, new_data1, n_rows1)
    # (4) SAMPLER HDU'S HEADER AND DATA
    #vegas_fits = write_tbl(vegas_fits, 4, new_data4, n_rows4)
    #vegas_fits = write_hdr(vegas_fits, 4, new_hdr4)
    # (6) DATA HDU'S DATA
    vegas_fits = write_tbl(vegas_fits, 6, new_data6, n_rows6)
    # SAVE THE NEW FILE
    write_file(vegas_fits, new_path, overwrite_bool=True)
    # CLOSE THE ORIGINAL FILE
    vegas_fits.close()

def write_new_SDFITS(session):
    fpath_original = f'/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/original/{session}/SDFITS/'
    fpath_new = f'/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/{session}/SDFITS/'
    # OPEN THE ORIGINAL FILES
    for spath in glob.glob(fpath_original + '*.fits'):
        old_sdfits = open_original_fits(old_path)
    # (0) PRIMARY HDU'S HEADER
    old_sdfits = write_hdr(old_sdfits, 0, new_hdr0)
    # (1) SINGLE DISH HDU'S DATA
    old_data1 = old_sdfits[1].data['DATA']
    old_indx1 = old_sdfits[1].data['SCAN'] == scannum
    old_data1[old_indx1] = new_data1
    old_sdfits = write_tbl(old_sdfits, 1, old_data1, n_rows1)
    old_sdfits = write_hdr(old_sdfits, 1, new_hdr1)
    # SAVE THE NEW FILE
    write_file(old_sdfits, new_path, overwrite_bool=True)
    # CLOSE THE ORIGINAL FILE
    old_sdfits.close()