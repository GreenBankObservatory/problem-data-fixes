import pandas as pd
import numpy as np
import glob
from scipy.ndimage import shift
from tqdm import tqdm
from core import *
from load_session import *
from write_fits import *

def fix_lo1a_vegas(datapath, session_code):
    lo1a_fnames = glob.glob(f"{datapath}/{session_code}/LO1A/*.fits")
    banks = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    for lpath_original in tqdm(lo1a_fnames, leave=False):
        lpath_short = lpath_original.split('/')[-1].replace(".fits", "")
        n_banks = len(glob.glob(f"{datapath}/{session_code}/VEGAS/{lpath_short}*.fits"))
        for bi in tqdm(range(n_banks)):
            bank = banks[bi]
            lpath_modified = lpath_original.replace("original", "modified")
            vpath_modified = lpath_modified.replace("LO1A", "VEGAS").replace(".fits", f"{bank}.fits")
            vpath_original = vpath_modified.replace("modified", "original")
            session = GBTSession(session_code, vpath_original)

            vdata = session.VEGAS_DATA
            vdata_shape = np.shape(vdata)

            vegas_new_hdr0 = {'PD_NAME':'23_EPHEM'}
            #vegas_new_data1 = 
            #vegas_n_rows1 = 
            #vegas_new_hdr4 =  
            #vegas_new_data4 = 
            #vegas_n_rows4 = 
            vegas_new_tbl6 = np.empty(np.shape(session.VEGAS_DATA)) # Needs filling
            vegas_n_rows6 = len(session.DMJD) # Good to go

            lo1a_new_hdr0 = {'PD_NAME':'23_EPHEM'}
            lo1a_new_data3 = {
                'DMJD':np.empty(vegas_n_rows6),
                'RA':np.empty(vegas_n_rows6),
                'DEC':np.empty(vegas_n_rows6),
                'LO1FREQ':np.empty(vegas_n_rows6),
                'VFRAME':np.empty(vegas_n_rows6),
                'RVSYS':np.empty(vegas_n_rows6)
            }
            lo1a_n_rows3 = len(session.DMJD)

            for i in tqdm(range(vegas_n_rows6), leave=False):
                # Calculate the correct VFRAME values
                vframe_i = calc_vframe(session.DMJD[i], session.RA[i], session.DEC[i], session.GO_VELDEF)

                # Calculate the actual VFRAME offset to apply 
                # given the offset that was originally applied
                vframe_off_i = calc_vframe_offset(session.LO1_DMJD, session.DMJD[i], session.LO1_VFRAME, vframe_i)

                # Calculate the frequency offset to shift all the data by
                f_offset = calc_f_offset(session.LO1_RESTFRQ, vframe_off_i)
                
                # Calculate the number of channels to shift by
                channel_shift_i = calc_channel_offset(f_offset, session.VEGAS_CHANBW[0])
                
                # Calculate the LO1FREQ
                lo1freq_i = sky2lo(session.LO1_RESTFRQ, session.LO1_LOMULT, session.LO1_IFFREQ, vframe_i, session.LO1_SIDEBAND, 'rel', session.LO1_S_VEL)

                # Calculate the RVSYS
                rvsys_i = calc_rvsys(session.LO1_S_VEL, vframe_i)

                # Add integration's shifted data to the new VEGAS table
                vegas_new_tbl6[i] = shift(session.VEGAS_DATA[i], [0,0,channel_shift_i], cval=np.NaN)
                # Add integration's information to the LO1A table
                lo1a_new_data3['DMJD'][i] = session.DMJD[i]
                lo1a_new_data3['RA'][i] = session.RA[i]
                lo1a_new_data3['DEC'][i] = session.DEC[i]
                lo1a_new_data3['LO1FREQ'][i] = lo1freq_i
                lo1a_new_data3['VFRAME'][i] = vframe_i
                lo1a_new_data3['RVSYS'][i] = rvsys_i

            # Write a new VEGAS FITS file
            vegas_new_data6 = {'DATA':vegas_new_tbl6}
            write_new_VEGAS(session.VEGAS_FNAME, vpath_modified, vegas_new_hdr0, vegas_new_data6, vegas_n_rows6)
            # new_hdr0, new_data1, n_rows1, new_hdr4, new_data4, n_rows4, 

            # Write a new LO1A FITS file
            write_new_LO1A(lpath_original, lpath_modified, lo1a_new_hdr0, lo1a_new_data3, lo1a_n_rows3)