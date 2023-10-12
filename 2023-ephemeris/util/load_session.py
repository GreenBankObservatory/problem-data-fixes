import os
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
import pandas as pd

class GBTSession:
    def __init__(self, LOADPATH, SESSION, LO1A_NAME, VEGAS_NAME):
        self.LOADPATH = LOADPATH
        self.SESSION = SESSION
        self.LO1A_NAME = LO1A_NAME
        self.VEGAS_NAME = VEGAS_NAME
        self.load_VEGAS()
        self.load_LO1()
        self.load_Antenna()
        self.load_GO()

    def load_VEGAS(self):
        ''' Loads the necessary information from the VEGAS FITS files '''

        # OPEN FILES
        self.VEGAS_PATH = os.path.join(self.LOADPATH, self.SESSION, "VEGAS", self.VEGAS_NAME)
        data = fits.open(self.VEGAS_PATH)

        # HDU 0 (PRIMARY)
        self.VEGAS_SCANNUM = data[0].header['SCAN']
        self.VEGAS_BANK = data[0].header['BANK'].strip()         # Bank (A-H)
        self.VEGAS_NCHAN = int(data[0].header['NCHAN'])          # Number of channels
        self.VEGAS_NPOL = int(data[0].header['NPOL'])            # Number of polarizations
        self.VEGAS_NSUBBAND = int(data[0].header['NSUBBAND'])    # Number of subbands (windows)
        self.VEGAS_ADCSAMPF = int(data[0].header['ADCSAMPF'])    # ADC Sampling frequency in Hz
        self.VEGAS_SUBFREQ = {}                                  # Center frequencies of each subband
        for i in range(self.VEGAS_NSUBBAND):
            self.VEGAS_SUBFREQ[f'SUB{i}FREQ'] = float(data[0].header[f'SUB{i}FREQ'])
        self.VEGAS_BW = data[0].header['BASE_BW']                # Overall bandwidth (MHz)

        # HDU 6 (DATA)
        self.DMJD = data[6].data['DMJD']                   # DMJD of integration   | dim: (n, )
        self.VEGAS_DATA = data[6].data['DATA']             # Counts                | dim: (n, SAMPLER, ACT_STATE, n_chan)
        self.VEGAS_NROWS = len(self.DMJD)

        # HDU 1 (SPURS)
        self.VEGAS_SPURCHAN = data[1].data['SPURCHAN']           # Channel index of spur
        self.VEGAS_SPURFREQ = data[1].data['SPURFREQ']           # Frequency of spur (Hz)

        # HDU 4 (SAMPLER)
        self.VEGAS_CRPIX1 = data[4].header['CRPIX1']             # Index of reference freq's channel
        self.VEGAS_CRVAL1 = data[4].data['CRVAL1']               # Channel center frequencies
        self.VEGAS_CDELT1 = pd.Series([data[4].data['CDELT1'] for i in range(len(self.DMJD)) ])       # Channel bandwidths

        # CLOSE FILE
        data.close()

    def load_LO1(self):
        ''' Loads the necessary information from the LO1A FITS files '''
        
        # OPEN FILES
        self.LO1A_PATH = os.path.join(self.LOADPATH, self.SESSION, "LO1A", self.LO1A_NAME)
        data = fits.open(self.LO1A_PATH)

        # HDU 1 (PRIMARY)
        self.LO1_RESTFRQ    = pd.Series([data[0].header['RESTFRQ'] for i in range(len(self.DMJD)) ])
        self.LO1_IFFREQ     = pd.Series([data[0].header['IFFREQ'] for i in range(len(self.DMJD)) ])
        self.LO1_LOOFFSET   = pd.Series([data[0].header['LOOFFSET'] for i in range(len(self.DMJD)) ])
        self.LO1_LOMULT     = pd.Series([data[0].header['LOMULT'] for i in range(len(self.DMJD)) ])
        self.LO1_REQDPTOL   = data[0].header['REQDPTOL']
        self.LO1_SIDEBAND   = pd.Series([data[0].header['SIDEBAND'] for i in range(len(self.DMJD)) ])

        # HDU 2 (SOUVEL)
        self.LO1_S_VEL      = pd.Series([data[2].data['VELOCITY'] for i in range(len(self.DMJD)) ])  # m/s

        # HDU 3 (LO1TBL)
        self.LO1_DMJD    = data[3].data['DMJD']
        self.LO1_RA      = self.interp_LO1A(data[3].data['RA'])
        self.LO1_DEC     = self.interp_LO1A(data[3].data['DEC'])
        self.LO1_LO1FREQ = self.interp_LO1A(data[3].data['LO1FREQ'])
        self.LO1_VFRAME  = self.interp_LO1A(data[3].data['VFRAME'])
        self.LO1_RVSYS   = self.interp_LO1A(data[3].data['RVSYS'])
        self.LO1_DMJD    = self.interp_LO1A(self.LO1_DMJD)

        # CLOSE FILE
        data.close()

    def interp_LO1A(self, l_row):
        """ Interpolare the LO1 table to have the same number of rows as VEGAS data """
        if len(l_row) == 1:
            l_row_interp = np.ones(self.VEGAS_NROWS)*l_row[0]
        else:
            f_l_row = interp1d(self.LO1_DMJD, l_row, kind='nearest')
            l_row_interp = f_l_row(self.DMJD)
        return pd.Series(l_row_interp)

    def load_Antenna(self):
        ''' Loads the necessary information from the Antenna FITS files '''

        # OPEN FILES
        self.ANTENNA_PATH = f"/home/gbtdata/{self.SESSION}/Antenna/{self.LO1A_NAME}"
        data = fits.open(self.ANTENNA_PATH)

        # HDU 2 (ANTPOSPF)
        ANT_DMJD    = data[2].data['DMJD']      # DMJD
        ANT_RA      = data[2].data['RAJ2000']   # RA (deg J2000)
        ANT_DEC     = data[2].data['DECJ2000']  # Dec (deg J2000)

        f_RA = interp1d(ANT_DMJD, ANT_RA, kind='linear')
        f_DEC = interp1d(ANT_DMJD, ANT_DEC, kind='linear')
        self.RA = pd.Series(f_RA(self.DMJD))
        self.DEC = pd.Series(f_DEC(self.DMJD))

        # CLOSE FILE
        data.close()

    def load_GO(self):
        ''' Loads the necessary information from the GO FITS files '''

        # OPEN FILES
        self.GO_PATH = f"/home/gbtdata/{self.SESSION}/GO/{self.LO1A_NAME}"
        data = fits.open(self.GO_PATH)

        # HDU 0 (PRIMARY)
        self.GO_VELDEF = pd.Series([data[0].header['VELDEF'] for i in range(len(self.DMJD))])
        
        # CLOSE FILE
        data.close()