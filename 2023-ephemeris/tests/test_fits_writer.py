from core.core import *
from load_session import *
from write_fits import *
import os

def test_LO1A_writer():

    original_fits = open_original_fits("./tests/local-data/sample_LO1A.fits")
    n_hdu = get_n_hdu(original_fits)

    hdr0_dict = {     
        'OBJECT': 'FITS TEST',                       
        'ADDME': 'ADDED!',
    }

    data3_dict = {
        'DMJD': [0, 1, 2, 3],
        'RA': [0, 1, 2, 3],
        'DEC': [0, 1, 2, 3],
        'LO1FREQ': [0, 1, 2, 3],
        'VFRAME': [0, 1, 2, 3],
        'RVSYS': [0, 1, 2, 3]
    }

    original_fits = write_hdr(original_fits, 0, hdr0_dict)
    original_fits = write_tbl(original_fits, 3, data3_dict, n_rows=4)
    write_file(original_fits, './tests/local-data/new_LO1A.fits', overwrite_bool=True)
    original_fits.close()

    new_fits = fits.open('./tests/local-data/new_LO1A.fits')
    
    # Make sure everything got added
    assert new_fits[0].header['OBJECT'] == 'FITS TEST'
    assert new_fits[0].header['ADDME'] == 'ADDED!'
    assert new_fits[3].header['EXTNAME'] == 'LO1TBL'
    assert new_fits[3].data['DMJD'][3] == 3
    new_fits.close()
    os.remove('./tests/local-data/new_LO1A.fits')

    # Make sure the original didn't get overwritten
    old_fits = fits.open('./tests/local-data/sample_LO1A.fits')
    assert old_fits[0].header['OBJECT'] == 'TEST'
    old_fits.close()