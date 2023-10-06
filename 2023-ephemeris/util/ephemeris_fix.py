import pandas as pd
import numpy as np
import os, glob
from scipy.ndimage import shift
from tqdm import tqdm
from astropy.time import Time

# LOCAL IMPORTS
from .core import *
from .tools import *
from .load_session import *
from .fits import *

#from rich.progress import TaskID

def get_lo1a_files(loadpath, session):
    """ Get the names of all LO1A files for the session """
    l_paths = list(glob.glob(f"{loadpath}/{session}/LO1A/*.fits"))
    l_names = [os.path.basename(lp) for lp in l_paths]
    return l_paths, l_names

def get_vegas_files(loadpath, session, l_name):
    """ Find all VEGAS files associated with an LO1A file """
    v_names_search = l_name.replace(".fits", "*.fits")
    v_paths = list(glob.glob(f"{loadpath}/{session}/VEGAS/{v_names_search}"))
    v_names = [os.path.basename(vp) for vp in v_paths]
    return v_paths, v_names

def make_lo1a_dict(session):
    """ Make a dictionary to hold the new LO1A data """
    n_rows = session.VEGAS_NROWS
    lo1a_dict = {
        0:{
            'hdr':{
                'HISTORY': 'File modified as part of Modification Request 5Q323',   # [12.2]
                'HISTORY': 'Original data can be found in problem-data/',           # [12.2]
                'REQDPTOL': str(session.VEGAS_FREQRES[0])                                   # [12.1]
                }
            }, 
        3:{
            'data':{
                'DMJD':np.empty(n_rows),
                'RA':np.empty(n_rows),
                'DEC':np.empty(n_rows),
                'LO1FREQ':np.empty(n_rows),
                'VFRAME':np.empty(n_rows),
                'RVSYS':np.empty(n_rows)
                }
            }
        }
    return lo1a_dict

def make_vegas_dict(session):
    """ Make a dictionary to hold the new VEGAS data """
    vegas_dict = {
        0:{
            'hdr':{
                'HISTORY': 'File modified as part of Modification Request 5Q323',   # [12.2]
                'HISTORY': 'Original data can be found in problem-data/'            # [12.2]
                }
            }, 
        1:{
            'hdr':{
                'COMMENT': 'The known spurs are calculated from the formula',
                'COMMENT': 'SPUR_CHAN = (J*ADCSAMPF/64-CRVAL1)/CDELT1+CRPIX1',
                'COMMENT': 'Where J ranges from 0 to 32',
                'COMMENT': 'Spurs outside the bandpass are excluded from the table',
                'COMMENT': 'Note the center channel (CRPIX1) is always flagged as a spur'
                },
            'data':{
                'SPURCHAN':session.VEGAS_SPURCHAN,
                'SPURFREQ':session.VEGAS_SPURFREQ
                }
            }, 
        4:{
            'data':{
                'CRVAL1':session.VEGAS_CRVAL1,
                }
            }, 
        6:{
            'data':{'data':session.VEGAS_DATA}
            }
        }
    return vegas_dict

def calc_progress(loadpath, session_name):
    l_paths, l_names = get_lo1a_files(loadpath, session_name)
    n_progress = len(l_names) + 1
    for l_name in l_names:
        v_paths, v_names = get_vegas_files(loadpath, session_name, l_name)
        n_progress += len(v_names)*4
    return n_progress

def main_fix(progress, task_id, loadpath, savepath, session_name):
    """ Runs the ephemeris fix """
    # Estimate amount of time
    progress.start_task(task_id)
    n_progress = calc_progress(loadpath, session_name)
    progress.update(task_id, total=n_progress)
    progress.update(task_id, advance=1)

    # [1] Create a GBT object
    gbt = create_obj_gbt()

    # Iterate through LO1A files
    l_paths, l_names = get_lo1a_files(loadpath, session_name)
    for k, l_name in enumerate(l_names):
        l_path_og = l_paths[k]
        l_path_save = os.path.join(savepath, session_name, "LO1A", l_name)
        v_paths, v_names = get_vegas_files(loadpath, session_name, l_name)
        for j, v_name in enumerate(v_names):
            v_path_og = v_paths[j]
            v_path_save = os.path.join(savepath, session_name, "VEGAS", v_name)
            
            # [2] Load the session
            session = GBTSession(loadpath, session_name, l_name, v_name)
            n_rows = len(session.DMJD)
            lo1a_dict = make_lo1a_dict(session)
            vegas_dict = make_vegas_dict(session)
            progress.update(task_id, advance=1)
            for i in range(n_rows):
                # [3] Calculate the VFRAME (m/s) for each VEGAS data row
                vframe_i = calc_vframe(gbt, session.DMJD[i], session.RA[i], session.DEC[i], session.GO_VELDEF)
                lo1a_dict[3]['data']['VFRAME'][i] = vframe_i # [4.7]
                # [4] Calculate the VFRAME_OFFSET
                vframe_off_i = calc_vframe_offset(session.LO1_VFRAME[i], vframe_i)
                # [5] Calculate a frequency offset (F_OFFSET) from VFRAME_OFFSET
                f_offset_i = calc_f_offset(session.LO1_RESTFRQ, vframe_off_i)
                # [6] Calculate a channel shift from the frequency offset
                channel_shift_i = calc_channel_offset(f_offset_i, session.VEGAS_CDELT1[0], session.LO1_SIDEBAND)
                # [7] Calculate the new LO1 frequencies
                lo1freq_i = sky2lo(session.LO1_RESTFRQ, session.LO1_LOMULT, session.LO1_IFFREQ, vframe_i, session.LO1_SIDEBAND, 'rel', session.LO1_S_VEL)
                lo1a_dict[3]['data']['LO1FREQ'][i] = lo1freq_i # [7.3]
                # [8] Calculate the new RVSYS values
                rvsys_i = calc_rvsys(session.LO1_S_VEL, vframe_i, session.GO_VELDEF)
                lo1a_dict[3]['data']['RVSYS'][i] = rvsys_i # [8.3]
                # [9] Shift the rows of data in the VEGAS FITS files by their channel offsets
                vegas_dict[6]['data']['data'][i] = shift(vegas_dict[6]['data']['data'][i], [0,0,channel_shift_i], cval=np.NaN)
                # [10] Write the DMJD, RA, and DEC values in the row of the LO1TBL
                lo1a_dict[3]['data']['DMJD'][i] = session.DMJD[i]
                lo1a_dict[3]['data']['RA'][i] = session.RA[i]
                lo1a_dict[3]['data']['DEC'][i] = session.DEC[i]
            
            progress.update(task_id, advance=1)
                
            # [11] Shift the VEGAS spur locations
            n_cri = np.shape(session.VEGAS_CRVAL1)[0]
            for cri in range(n_cri):
                vegas_dict[4]['data']['CRVAL1'][cri] = float(vegas_dict[4]['data']['CRVAL1'][cri] + f_offset_i)
                vegas_dict[1]['data']['SPURCHAN'][cri] = int((cri*session.VEGAS_ADCSAMPF/64-vegas_dict[4]['data']['CRVAL1'][cri])/session.VEGAS_CDELT1[cri]+session.VEGAS_CRPIX1)
                vegas_dict[1]['data']['SPURFREQ'][cri] = float(vegas_dict[1]['data']['SPURFREQ'][cri]+f_offset_i)
            progress.update(task_id, advance=1)

            # Writing the new VEGAS files
            write_new_VEGAS(v_path_og, v_path_save, vegas_dict)
            progress.update(task_id, advance=1)

        # Writing the new LO1A files
        write_new_LO1A(l_path_og, l_path_save, lo1a_dict)
        progress.update(task_id, advance=1)
            
