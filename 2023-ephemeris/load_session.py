import glob
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.interpolate import interp1d

class GBTSession:
    def __init__(self, SESSION, VEGAS_FNAME):
        self.SESSION = SESSION
        self.VEGAS_FNAME = VEGAS_FNAME
        self.load_VEGAS()
        self.load_LO1()
        self.load_Antenna()
        self.load_GO()

    def load_VEGAS(self):
        ''' Loads the necessary information from the VEGAS FITS files '''

        # OPEN FILES
        self.FNAME_SHORT = self.VEGAS_FNAME.split('/')[-1][:-6]
        data = fits.open(self.VEGAS_FNAME)

        # HDU 0 (PRIMARY)
        self.VEGAS_SCANNUM = data[0].header['SCAN']
        self.VEGAS_BANK = data[0].header['BANK'].strip()         # Bank (A-H)
        self.VEGAS_NCHAN = int(data[0].header['NCHAN'])          # Number of channels
        self.VEGAS_NPOL = int(data[0].header['NPOL'])            # Number of polarizations
        self.VEGAS_NSUBBAND = int(data[0].header['NSUBBAND'])    # Number of subbands (windows)
        self.VEGAS_SUBFREQ = []                                  # Center frequencies of each subband
        for i in range(self.VEGAS_NSUBBAND):
            self.VEGAS_SUBFREQ.append(float(data[0].header[f'SUB{i}FREQ']))
        self.VEGAS_BW = data[0].header['BASE_BW']                # Overall bandwidth (MHz)

        # HDU 1 (SPURS)
        self.VEGAS_SPURCHAN = data[1].data['SPURCHAN']           # Channel index of spur
        self.VEGAS_SPURFREQ = data[1].data['SPURFREQ']           # Frequency of spur (Hz)

        # HDU 4 (SAMPLER)
        self.VEGAS_CRPIX = data[4].header['CRPIX1']              # Index of reference freq's channel
        self.VEGAS_CHANFREQ = data[4].data['CRVAL1']             # Channel center frequencies
        self.VEGAS_CHANBW = data[4].data['CDELT1']               # Channel bandwidths

        # HDU 6 (DATA)
        self.DMJD = data[6].data['DMJD']                   # DMJD of integration   | dim: (n, )
        self.VEGAS_DATA = data[6].data['DATA']             # Counts                | dim: (n, SAMPLER, ACT_STATE, n_chan)

        # CLOSE FILE
        data.close()

    def load_LO1(self):
        ''' Loads the necessary information from the LO1A FITS files '''
        
        # OPEN FILES
        FNAME_LO1 = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/original/{self.SESSION}/LO1A/{self.FNAME_SHORT}.fits"
        data = fits.open(FNAME_LO1)

        # HDU 1 (PRIMARY)
        self.LO1_RESTFRQ    = data[0].header['RESTFRQ']
        self.LO1_IFFREQ     = data[0].header['IFFREQ']
        self.LO1_LOOFFSET   = data[0].header['LOOFFSET']
        self.LO1_LOMULT     = data[0].header['LOMULT']
        self.LO1_REQDPTOL   = data[0].header['REQDPTOL']
        self.LO1_SIDEBAND   = data[0].header['SIDEBAND']

        # HDU 2 (SOUVEL)
        self.LO1_S_DMJD     = data[2].data['DMJD']
        self.LO1_S_VEL      = data[2].data['VELOCITY']  # m/s
        self.LO1_S_VDOT     = data[2].data['VDOT']      # m/s/s
        self.LO1_S_VDOTDOT  = data[2].data['VDOTDOT']   # m/s/s/s

        # HDU 3 (LO1TBL)
        self.LO1_DMJD    = data[3].data['DMJD']
        self.LO1_RA      = data[3].data['RA']
        self.LO1_DEC     = data[3].data['DEC']
        self.LO1_LO1FREQ = data[3].data['LO1FREQ']
        self.LO1_VFRAME  = data[3].data['VFRAME']
        self.LO1_RVSYS   = data[3].data['RVSYS']

        # CLOSE FILE
        data.close()

    def load_Antenna(self):
        ''' Loads the necessary information from the Antenna FITS files '''

        # OPEN FILES
        FNAME_ANTENNA = f"/home/gbtdata/{self.SESSION}/Antenna/{self.FNAME_SHORT}.fits"
        data = fits.open(FNAME_ANTENNA)

        # HDU 2 (ANTPOSPF)
        ANT_DMJD    = data[2].data['DMJD']      # DMJD
        ANT_RA      = data[2].data['RAJ2000']   # RA (deg J2000)
        ANT_DEC     = data[2].data['DECJ2000']  # Dec (deg J2000)

        f_RA = interp1d(ANT_DMJD, ANT_RA)
        f_DEC = interp1d(ANT_DMJD, ANT_DEC)
        self.RA = f_RA(self.DMJD)
        self.DEC = f_DEC(self.DMJD)

        # CLOSE FILE
        data.close()

    def load_GO(self):
        ''' Loads the necessary information from the GO FITS files '''

        # OPEN FILES
        FNAME_GO = f"/home/gbtdata/{self.SESSION}/GO/{self.FNAME_SHORT}.fits"
        data = fits.open(FNAME_GO)

        # HDU 0 (PRIMARY)
        self.GO_VELDEF = data[0].header['VELDEF']

        # CLOSE FILE
        data.close()