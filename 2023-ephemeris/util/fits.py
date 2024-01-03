from astropy.io import fits
import numpy as np
import os, datetime

def open_original_fits(fits_path):
    """ Open the original FITS file """
    original_fits = fits.open(fits_path)
    return original_fits

def check_hdr_key(k, hdu_num):
    """ Check if a key should be added to the new header """
    keep_key = True
    keys_to_ignore = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 
                      'NAXIS2', 'PCOUNT', 'GCOUNT', 'TFIELDS']
    starts_to_ignore = ['COMMENT']#['TTYPE', 'TFORM', 'TUNIT', 'COMMENT']
    if k in keys_to_ignore:
        keep_key = True
    for s in starts_to_ignore:
        if k.startswith(s):
            keep_key = False
    return keep_key
    
def get_tbl_hdr(tbl, hdu_num):
    """ Get the header of a given BinTableHDU """
    all_keys = {}
    for k in tbl.header.keys():
        if check_hdr_key(k, hdu_num):
            all_keys[k] = {
                'val': tbl.header[k],
                'comment': tbl.header.comments[k]
            }
    return all_keys

def write_tbl(data_hdulist, hdu_num, new_data, n_rows):
    """ Modify a BinTableHDU's data with new values """
    old_tbl = data_hdulist[hdu_num]
    new_colnames = new_data.keys()
    new_tbl = fits.BinTableHDU.from_columns(old_tbl.columns, nrows=n_rows)
    new_hdr = get_tbl_hdr(old_tbl, hdu_num)
    for new_col in new_colnames:
        new_tbl.data[new_col] = new_data[new_col]
        #new_tbl.data.set(new_col, new_data[new_col]['val'], new_data[new_col]['comment'])
    for new_key in new_hdr.keys():
        #new_tbl.header.set(new_key, new_hdr[new_key])
        new_tbl.header.set(new_key, new_hdr[new_key]['val'], new_hdr[new_key]['comment'])
    data_hdulist[hdu_num] = new_tbl
    return data_hdulist

def write_hdr(data_hdulist, hdu_num, new_hdr):
    """ Modify an HDU's header with new values """
    for k in new_hdr.keys():
        new_k = str(new_hdr[k])
        data_hdulist[hdu_num].header.set(k, new_k)
    for i in range(8):
        try:
            k = f"SUB{i}FREQ"
            new_k = str(int(new_hdr[k])) + '.'
            data_hdulist[hdu_num].header.set(k, new_k)
        except:
            pass
    return data_hdulist

def write_file(lock, fits_file, save_path):
    """ Write the modified FITS file to a new location """

    t_now = datetime.datetime.now()
    fits_file[0].header['HISTORY'] = 'File modified as part of Modification Request 5Q323'
    fits_file[0].header['HISTORY'] = 'Original data can be found in problem-data/'
    fits_file[0].header['HISTORY'] = f'Last modified: {t_now.strftime("%Y-%m-%d %H:%M:%S")}'
    fits_file.writeto(save_path, overwrite=True, output_verify='ignore')
    '''
    except:
        with lock:
            print(f"ERROR SAVING {fits_file} TO {save_path}.")
            print(os.path.isfile(save_path), save_path)
            print(os.path.isdir(os.path.dirname(save_path)), os.path.dirname(save_path))
            print(os.path.isdir(os.path.dirname(os.path.dirname(save_path))), os.path.dirname(os.path.dirname(save_path)))
            print(os.path.isdir(os.path.dirname(os.path.dirname(os.path.dirname(save_path)))), os.path.dirname(os.path.dirname(os.path.dirname(save_path))))
            fits_file.writeto("/home/sandboxes/vcatlett/ERRORFILE.fits", overwrite=True, output_verify='ignore')
        '''

def write_new_LO1A(lock, old_path, new_path, lo1a_dict):
    """ Write a new LO1A FITS file """
    lo1a_fits = open_original_fits(old_path)
    lo1a_fits = write_hdr(lo1a_fits, 0, lo1a_dict[0]['hdr'])
    n_rows3 = len(lo1a_dict[3]['data']['DMJD'])
    lo1a_fits = write_tbl(lo1a_fits, 3, lo1a_dict[3]['data'], n_rows3)
    write_file(lock, lo1a_fits, new_path)
    lo1a_fits.close()

def write_new_VEGAS(lock, old_path, new_path, vegas_dict):
    """ Write a new VEGAS FITS file """
    vegas_fits = open_original_fits(old_path)
    vegas_fits = write_hdr(vegas_fits, 0, vegas_dict[0]['hdr'])
    vegas_fits = write_hdr(vegas_fits, 1, vegas_dict[1]['hdr'])
    n_rows1 = len(vegas_dict[1]['data']['SPURFREQ'])
    vegas_fits = write_tbl(vegas_fits, 1, vegas_dict[1]['data'], n_rows1)
    n_rows6 = np.shape(vegas_dict[6]['data']['data'])[0]
    vegas_fits = write_tbl(vegas_fits, 6, vegas_dict[6]['data'], n_rows6)
    vegas_fits[1].header['COMMENT'] = 'The known spurs are calculated from the formula'
    vegas_fits[1].header['COMMENT'] = 'SPUR_CHAN = (J*ADCSAMPF/64-CRVAL1)/CDELT1+CRPIX1'
    vegas_fits[1].header['COMMENT'] = 'Where J ranges from 0 to 32'
    vegas_fits[1].header['COMMENT'] = 'Spurs outside the bandpass are excluded from the table'
    vegas_fits[1].header['COMMENT'] = 'Note the center channel (CRPIX1) is always flagged as a spur'
    write_file(lock, vegas_fits, new_path)
    vegas_fits.close()

def fill_sdfits(loadpath, savepath, session_code):
    """ Fill a new SDFITS file """
    pass