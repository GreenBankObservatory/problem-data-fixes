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
        self.load_IF()
        self.load_LO1()
        self.load_GO()
        self.load_Antenna()

    def load_VEGAS(self):
        ''' Loads the necessary information from the VEGAS FITS files '''

        # OPEN FILES
        self.VEGAS_PATH = os.path.join(self.LOADPATH, self.SESSION, "VEGAS", self.VEGAS_NAME)
        data = fits.open(self.VEGAS_PATH, memmap=True, cache=False, lazy_load_hdus=True)

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
        self.DMJD = list(data[6].data['DMJD'])                   # DMJD of integration   | dim: (n, )
        self.VEGAS_DATA = data[6].data['DATA']             # Counts                | dim: (n, SAMPLER, ACT_STATE, n_chan)
        self.VEGAS_NROWS = len(self.DMJD)
        self.VEGAS_INTEGNUM = data[6].data['INTEGNUM']

        # HDU 1 (SPURS)
        self.VEGAS_SPURCHAN = data[1].data['SPURCHAN']           # Channel index of spur
        self.VEGAS_SPURFREQ = data[1].data['SPURFREQ']           # Frequency of spur (Hz)

        # HDU 4 (SAMPLER)
        self.VEGAS_CRPIX1 = data[4].header['CRPIX1']            # Index of reference freq's channel
        self.VEGAS_CRVAL1 = pd.Series([data[4].data['CRVAL1'][0] for i in range(len(self.DMJD)) ])               # Channel center frequencies
        self.VEGAS_CDELT1 = pd.Series([data[4].data['CDELT1'][0] for i in range(len(self.DMJD)) ])       # Channel bandwidths
        self.VEGAS_SIGREF = data[3].data['SIGREF']
        self.VEGAS_CAL = data[3].data['CAL']
        # CLOSE FILE
        data.close()
        del data

    def load_LO1(self):
        ''' Loads the necessary information from the LO1A FITS files '''
        
        # OPEN FILES
        self.LO1A_PATH = os.path.join(self.LOADPATH, self.SESSION, "LO1A", self.LO1A_NAME)
        data = fits.open(self.LO1A_PATH, memmap=True, cache=False, lazy_load_hdus=True)

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
        del data

    def load_IF(self):
        ''' Loads the necessary information from the LO1A FITS files '''
        
        # OPEN FILES
        self.IF_PATH = f"/home/gbtdata/{self.SESSION}/IF/{self.LO1A_NAME}"
        data = fits.open(self.IF_PATH, memmap=True, cache=False, lazy_load_hdus=True)

        self.IF_FEED = data[1].data["FEED"]
        self.IF_POL = data[1].data["POLARIZE"]

        # CLOSE FILE
        data.close()
        del data

    def interp_LO1A(self, l_row):
        """ Interpolate the LO1 table to have the same number of rows as VEGAS data """
        if len(l_row) == 1:
            l_row_interp = np.ones(self.VEGAS_NROWS)*l_row[0]
        else:
            f_l_row = interp1d(self.LO1_DMJD, l_row, kind='nearest', bounds_error=None, fill_value="extrapolate")
            l_row_interp = f_l_row(self.DMJD)
        return pd.Series(l_row_interp)

    def load_Antenna(self):
        ''' Loads the necessary information from the Antenna FITS files '''

        self.RA = self.LO1_RA
        self.DEC = self.LO1_DEC

        try:
            # OPEN FILES
            self.ANTENNA_PATH = f"/home/gbtdata/{self.SESSION}/Antenna/{self.LO1A_NAME}"
            data = fits.open(self.ANTENNA_PATH, memmap=True, cache=False, lazy_load_hdus=True)

            # HDU 2 (ANTPOSPF)
            ANT_DMJD    = data[2].data['DMJD']      # DMJD
            ANT_RA      = data[2].data['RAJ2000']   # RA (deg J2000)
            ANT_DEC     = data[2].data['DECJ2000']  # Dec (deg J2000)

            indx_off = self.GO_SWSTATE.str.startswith("PSWITCHOFF")
            indx_on = ~indx_off
            #indx_on = list(indx_on[indx_on].index)
            #print("INDX_ON: ", indx_on)
            
            f_RA = interp1d(ANT_DMJD, ANT_RA, kind='linear', bounds_error=None, fill_value="extrapolate")
            f_DEC = interp1d(ANT_DMJD, ANT_DEC, kind='linear', bounds_error=None, fill_value="extrapolate")
            self.RA[indx_on] = list(f_RA(np.array(self.DMJD)[indx_on]))
            self.DEC[indx_on]= list(f_DEC(np.array(self.DMJD)[indx_on]))

            # CLOSE FILE
            data.close()
            del data
        except:
            pass

    def load_GO(self):
        ''' Loads the necessary information from the GO FITS files '''

        # OPEN FILES
        self.GO_PATH = f"/home/gbtdata/{self.SESSION}/GO/{self.LO1A_NAME}"
        data = fits.open(self.GO_PATH, memmap=True, cache=False, lazy_load_hdus=True)

        # HDU 0 (PRIMARY)
        self.GO_VELDEF = pd.Series([data[0].header['VELDEF'] for i in range(len(self.DMJD))])
        self.GO_SWSTATE = pd.Series([data[0].header['SWSTATE'] for i in range(len(self.DMJD))])
        
        # CLOSE FILE
        data.close()
        del data