from util.core import *
from util.load_session import *
from util.fits import *
import os

def test_fits_writer():
    """ Check that a FITS file can be opened, edited, and saved to a new file without changes to the original file """

    lock = None

    original_fits = open_original_fits("/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/pytest_data/sample_LO1A.fits")

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
    write_file(lock, original_fits, '/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/pytest_data/new_LO1A.fits')
    original_fits.close()

    new_fits = fits.open('/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/pytest_data/new_LO1A.fits')
    
    # Make sure everything got added
    assert new_fits[0].header['OBJECT'] == 'FITS TEST'
    assert new_fits[0].header['ADDME'] == 'ADDED!'
    assert new_fits[3].header['EXTNAME'] == 'LO1TBL'
    assert new_fits[3].data['DMJD'][3] == 3
    new_fits.close()
    os.remove('/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/pytest_data/new_LO1A.fits')

    # Make sure the original didn't get overwritten
    old_fits = fits.open('/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/pytest_data/sample_LO1A.fits')
    assert old_fits[0].header['OBJECT'] == 'TEST'
    old_fits.close()