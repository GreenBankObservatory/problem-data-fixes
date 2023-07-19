import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def calc_freqs(crpix1, cdelt1, crval1, nchan):
    indxs = np.array(range(1, nchan+1))
    freqs = np.add(np.multiply(np.subtract(indxs, crpix1), cdelt1), crval1)
    return freqs

def calc_tsys(tcal, calon, caloff):
    indx_low = np.round(len(calon)*0.1)
    indx_high = np.round(len(calon)*0.9)
    calon_80 = np.mean(calon[indx_low:indx_high])
    calonoff_80 = np.mean(np.subtract(calon[indx_low:indx_high], caloff[indx_low:indx_high]))
    tsys = np.add(tcal * (calon_80/calonoff_80), tcal/2)
    return tsys

def calc_tant(ref_tsys, sig_calon, sig_caloff, ref_calon, ref_caloff):
    sig = 0.5 * np.add(sig_calon, sig_caloff)
    ref = 0.5 * np.add(ref_calon, ref_caloff)
    tant = np.mult(ref_tsys, np.divide(np.subtract(sig, ref), ref))
    return tant

bad_fpath = "/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/original/AGBT21B_316_21/SDFITS/AGBT21B_316_21.raw.vegas.A.fits"
fixed_fpath = "/home/sandboxes/vcatlett/repos/github/problem-data-fixes/2023-ephemeris/tests/data/newVEGAS.fits"
good_fpath = "/home/sdfits/AGBT21B_316_26/AGBT21B_316_26.raw.vegas/AGBT21B_316_26.raw.vegas.A.fits"

bad_data = fits.open(bad_fpath)
fixed_data = fits.open(fixed_fpath)
good_data = fits.open(good_fpath)

bad_spectra = bad_data[1].data['DATA']
fixed_spectra = fixed_data[1].data['DATA']
good_spectra = good_data[1].data['DATA']

bad_scans = bad_data[1].data['SCAN']
fixed_scans = fixed_data[1].data['SCAN']
good_scans = good_data[1].data['SCAN']

bad_indx = bad_scans==100
fixed_indx = fixed_scans==100
good_indx = good_scans==6

bad_nchan = np.shape(bad_spectra)[1]
fixed_nchan = np.shape(fixed_spectra)[1]
good_nchan = np.shape(good_spectra)[1]

fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True)

cax = ax[0]
for i in range(len(bad_scans[bad_indx])):
    crpix1_i = bad_data[1].data['CRPIX1'][bad_indx][i]
    cdelt1_i = bad_data[1].data['CDELT1'][bad_indx][i]
    crval1_i = bad_data[1].data['CRVAL1'][bad_indx][i]
    bad_freqs = calc_freqs(crpix1_i, cdelt1_i, crval1_i, bad_nchan)
    bad_spec_i = bad_spectra[bad_indx][i]
    cax.plot(bad_freqs, bad_spec_i)

cax = ax[1]
for i in range(len(fixed_scans[fixed_indx])):
    crpix1_i = fixed_data[1].data['CRPIX1'][fixed_indx][i]
    cdelt1_i = fixed_data[1].data['CDELT1'][fixed_indx][i]
    crval1_i = fixed_data[1].data['CRVAL1'][fixed_indx][i]
    fixed_freqs = calc_freqs(crpix1_i, cdelt1_i, crval1_i, fixed_nchan)
    fixed_spec_i = fixed_spectra[fixed_indx][i]
    cax.plot(fixed_freqs, fixed_spec_i)

cax = ax[2]
for i in range(len(good_scans[good_indx])):
    crpix1_i = good_data[1].data['CRPIX1'][good_indx][i]
    cdelt1_i = good_data[1].data['CDELT1'][good_indx][i]
    crval1_i = good_data[1].data['CRVAL1'][good_indx][i]
    good_freqs = calc_freqs(crpix1_i, cdelt1_i, crval1_i, good_nchan)
    good_spec_i = good_spectra[good_indx][i]
    cax.plot(good_freqs, good_spec_i)

plt.show()

bad_data.close()
good_data.close()