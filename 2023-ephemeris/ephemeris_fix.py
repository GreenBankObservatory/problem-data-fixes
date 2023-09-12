import pandas as pd
import numpy as np
import glob
from scipy.ndimage import shift
from tqdm import tqdm
from core import *
from load_session import *
from write_fits import *
from astropy.time import Time

def closest(lst, K):
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return idx, lst[idx]


def fix_lo1a_vegas(datapath, session_code):
    lo1a_fnames = glob.glob(f"{datapath}/{session_code}/LO1A/*.fits")
    banks = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    dict_shifts = {
            'bank':[],
            'scan':[],
            'shift':[]
        }
    for lpath_original in lo1a_fnames:
        lpath_short = lpath_original.split('/')[-1].replace(".fits", "")
        n_banks = len(glob.glob(f"{datapath}/{session_code}/VEGAS/{lpath_short}*.fits"))
        for bi in range(n_banks):
            bank = banks[bi]
            lpath_modified = lpath_original.replace("original", "modified")
            vpath_modified = lpath_modified.replace("LO1A", "VEGAS").replace(".fits", f"{bank}.fits")
            vpath_original = vpath_modified.replace("modified", "original")
            #try:
            session = GBTSession(session_code, vpath_original)

            vdata = session.VEGAS_DATA
            vdata_shape = np.shape(vdata)

            #if len(vdata) > 0:
            #try:
            vegas_new_hdr0 = {'PD_NAME':'23_EPHEM'}
            spurchan = True
            try:
                vegas_new_data4 = {
                    'SPURCHAN':session.VEGAS_SPURCHAN
                }
                vegas_n_rows4 = np.shape(session.VEGAS_SPURCHAN)[0]
            except:
                spurchan = False
                print("Couldn't find SPURCHAN, but moving on. ")
            
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

            for i in range(vegas_n_rows6):
                # Calculate the correct VFRAME values
                vframe_i = calc_vframe(session.DMJD[i], session.RA[i], session.DEC[i], session.GO_VELDEF)

                # Calculate the actual VFRAME offset to apply 
                # given the offset that was originally applied
                vframe_off_i = calc_vframe_offset(session.LO1_DMJD, session.DMJD[i], session.LO1_VFRAME, vframe_i)

                # Calculate the frequency offset to shift all the data by
                f_offset = calc_f_offset(session.LO1_RESTFRQ, vframe_off_i)
                
                # Calculate the number of channels to shift by
                channel_shift_i = calc_channel_offset(f_offset, session.VEGAS_CHANBW[0], session.LO1_SIDEBAND)

                # Shift the spur channels by that much
                if spurchan: 
                    for si in range(np.shape(session.VEGAS_SPURCHAN)[0]):
                        vegas_new_data4['SPURCHAN'][si] += channel_shift_i 
                        if vegas_new_data4['SPURCHAN'][si] < 0:
                            vegas_new_data4['SPURCHAN'][si] = -999
                
                # Calculate the LO1FREQ
                lo1freq_i = sky2lo(session.LO1_RESTFRQ, session.LO1_LOMULT, session.LO1_IFFREQ, vframe_i, session.LO1_SIDEBAND, 'rel', session.LO1_S_VEL)

                # Calculate the RVSYS
                rvsys_i = calc_rvsys(session.LO1_S_VEL, vframe_i, session.GO_VELDEF)

                # Add integration's shifted data to the new VEGAS table
                vegas_new_tbl6[i] = shift(session.VEGAS_DATA[i], [0,0,channel_shift_i], cval=np.NaN)
                # Add integration's information to the LO1A table
                lo1a_new_data3['DMJD'][i] = session.DMJD[i]
                lo1a_new_data3['RA'][i] = session.RA[i]
                lo1a_new_data3['DEC'][i] = session.DEC[i]
                lo1a_new_data3['LO1FREQ'][i] = lo1freq_i
                lo1a_new_data3['VFRAME'][i] = vframe_i
                lo1a_new_data3['RVSYS'][i] = rvsys_i

                dict_shifts['bank'].append(bank)
                dict_shifts['scan'].append(session.VEGAS_SCANNUM)
                dict_shifts['shift'].append(channel_shift_i)
            # Write a new VEGAS FITS file
            vegas_new_data6 = {'DATA':vegas_new_tbl6}
            write_new_VEGAS(session.VEGAS_FNAME, vpath_modified, vegas_new_hdr0, vegas_new_data4, vegas_n_rows4,vegas_new_data6, vegas_n_rows6, spur_bool=spurchan)

            # Write a new LO1A FITS file
            write_new_LO1A(lpath_original, lpath_modified, lo1a_new_hdr0, lo1a_new_data3, lo1a_n_rows3)
            #except:
            #    print(f"Sorry, something went wrong with analyzing {vpath_modified}")
        #except:
            #    print(f"Sorry, something went wrong with opening {vpath_modified}")
    df_shifts = pd.DataFrame(dict_shifts)
    df_shifts.to_csv(f'{datapath.replace("original", "modified")}/{session_code}/row-shifts.csv', index=False)


def fix_sdfits(session_code):
    # PATHS TO THE FILES
    row_shift_path = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/{session_code}/row-shifts.csv"
    old_path = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/original/{session_code}/"
    new_path = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/{session_code}/"
    banks = ['A','B','C','D','E','F','G','H']
    # FIND ORIGINAL SDFITS FILES
    all_old_sdfits = glob.glob(old_path + 'SDFITS/*.fits')
   
    # OPEN THE INFO FILE
    row_shifts = pd.read_csv(row_shift_path)

    # COMPILE LO1A SCAN INFO
    all_lo1a = sorted(glob.glob(new_path + 'LO1A/*.fits'))
    all_vegas = sorted(glob.glob(new_path + 'VEGAS/*.fits'))
    lo1a_scans = {}
    for li in all_lo1a:
        data_li = fits.open(li)
        lo1a_scans[data_li[0].header['SCAN']] = li
        data_li.close()
    
    # COMPILE VEGAS SCAN INFO
    vegas_scans = {}
    for vi in all_vegas:
        data_vi = fits.open(vi)
        vbi = vi.split('.')[-2][-1]
        try:
            vegas_scans[int(data_vi[0].header['SCAN'])][vbi] = vi
        except:
            vegas_scans[int(data_vi[0].header['SCAN'])] = {vbi:vi}
        data_vi.close()

    # GO THROUGH SDFITS FILES
    for old_sd in all_old_sdfits:
        try:
            old_fits = fits.open(old_sd)
            old_data = old_fits[1].data['DATA']
            new_data = np.zeros(np.shape(old_data))
            #new_twarm = np.zeros(np.shape(old_fits[1].data['TWARM']))
            #new_tcold = np.zeros(np.shape(old_fits[1].data['TCOLD']))
            #new_tcal = np.zeros(np.shape(old_fits[1].data['TCAL']))
            new_vframes = []
            new_rvsyss = []
            new_hdr0 = {
                'PD_NAME':'2023-EPH'
            }
            old_scannums = old_fits[1].data['SCAN']
            old_dmjds = Time(old_fits[1].data['DATE-OBS']).mjd
            # GO THROUGH ALL THE SCANS
            old_sn = -10
            for si in range(len(old_scannums)):
                sn = old_scannums[si]
                # FIND THE RIGHT LO1A AND VEGAS FILES
                lo1a_fits = fits.open(lo1a_scans[sn])
                vegas_fits = fits.open(vegas_scans[sn][banks[old_fits[1].data['IFNUM'][si]]])
                
                lo1a_dmjds = lo1a_fits[3].data['DMJD']
                l_indx, close_lo1a = closest(lo1a_dmjds, old_dmjds[si])
                vegas_dmjds = vegas_fits[6].data['DMJD']
                v_indx, close_vegas = closest(vegas_dmjds, old_dmjds[si])

                n_shift = np.count_nonzero(np.isnan(vegas_fits[6].data['DATA'][v_indx][0][0]))
                sign_shift = np.isnan(vegas_fits[6].data['DATA'][v_indx][0][0][0]) # True = positive shift
                if not sign_shift:
                    n_shift = -n_shift
                
                new_row = shift(old_data[si], n_shift, cval=np.NaN)
                new_data[si] = new_row
                #new_twarm[si] = shift(new_twarm[si], n_shift, cval=np.NaN)
                #new_tcold[si] = shift(new_tcold[si], n_shift, cval=np.NaN)
                #new_tcal[si] = shift(new_tcal[si], n_shift, cval=np.NaN)

                new_vframe = lo1a_fits[3].data['VFRAME'][l_indx]

                new_vframes.append(new_vframe)
                new_rvsyss.append(lo1a_fits[3].data['RVSYS'][l_indx])

                lo1a_fits.close()
                vegas_fits.close()
            new_data_dict = {
                'DATA':new_data,
                'VFRAME':new_vframes,
                'RVSYS':new_rvsyss
                }
            write_new_SDFITS(old_sd, old_sd.replace('original', 'modified'), new_hdr0, new_data_dict, len(old_scannums))
        except:
            print(f"Something went wrong for {old_sd}")