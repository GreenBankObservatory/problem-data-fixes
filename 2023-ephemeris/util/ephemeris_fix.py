import pandas as pd
import numpy as np
import os, glob
from scipy.ndimage import shift
from tqdm import tqdm
from astropy.time import Time
import datetime

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

def make_lo1a_dict(session, t_now):
    """ Make a dictionary to hold the new LO1A data """
    n_rows = session.VEGAS_NROWS
    lo1a_dict = {
        0:{
            'hdr':{
                'HISTORY': 'File modified as part of Modification Request 5Q323',   # [12.2]
                'HISTORY': 'Original data can be found in problem-data/',           # [12.2]
                'HISTORY': f'Last modified: {t_now.strftime("%Y-%m-%d %H:%M:%S")}',
                'REQDPTOL': str(session.VEGAS_CDELT1[0])                            # [12.1]
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

def make_vegas_dict(session, t_now):
    """ Make a dictionary to hold the new VEGAS data """
    vegas_dict = {
        0:{
            'hdr':{
                'HISTORY': 'File modified as part of Modification Request 5Q323',   # [12.2]
                'HISTORY': 'Original data can be found in problem-data/',           # [12.2]
                'HISTORY': f'Last modified: {t_now.strftime("%Y-%m-%d %H:%M:%S")}'
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
                'CDELT1':session.VEGAS_CDELT1,
                }
            }, 
        6:{
            'data':{'data':session.VEGAS_DATA}
            }
        }
    for sf in session.VEGAS_SUBFREQ.keys():
        vegas_dict[0]['hdr'][sf] = session.VEGAS_SUBFREQ[sf]
    return vegas_dict

def calc_progress(loadpath, session_name):
    l_paths, l_names = get_lo1a_files(loadpath, session_name)
    n_progress = 1
    for l_name in l_names:
        v_paths, v_names = get_vegas_files(loadpath, session_name, l_name)
        if len(v_names) > 0:
            n_progress += len(v_names)*4 + 1
        else:
            n_progress += 1
    return n_progress

def check_empty(fpath):
    """ Check if a VEGAS FITS file is empty (no integrations) """
    file_empty = True
    fdata = fits.open(fpath)
    n_int = np.shape(fdata[6].data["DATA"])[0]
    if n_int > 0:
        file_empty = False
    return file_empty

def main_fix(lock, progress, overall_progress, task_id, loadpath, savepath, session_name):
    """ Runs the ephemeris fix """
    # Estimate amount of time
    progress.start_task(task_id)
    #n_progress = calc_progress(loadpath, session_name)
    #progress.update(task_id, total=n_progress)
    progress.update(task_id, advance=1)
    progress.update(overall_progress, advance=1)

    # [1] Create a GBT object
    gbt = create_obj_gbt(lock)
    t_now = datetime.datetime.now()

    # Iterate through LO1A files
    l_paths, l_names = get_lo1a_files(loadpath, session_name)
    for k, l_name in enumerate(l_names):
        l_path_og = l_paths[k]
        l_path_save = os.path.join(savepath, session_name, "LO1A", l_name)
        v_paths, v_names = get_vegas_files(loadpath, session_name, l_name)
        for j, v_name in enumerate(v_names):
            v_path_og = v_paths[j]
            v_path_save = os.path.join(savepath, session_name, "VEGAS", v_name)
            v_is_empty = check_empty(v_path_og)

            if v_is_empty:
                progress.update(task_id, advance=1)
                progress.update(overall_progress, advance=1)
            else:            
                # [2] Load the session
                session = GBTSession(loadpath, session_name, l_name, v_name)
                n_rows = len(session.DMJD)

                #with lock:
                #    print("\nSESSION: ", session_name, l_name, v_name)
                
                lo1a_dict = make_lo1a_dict(session, t_now)
                vegas_dict = make_vegas_dict(session, t_now)
                progress.update(task_id, advance=1)
                progress.update(overall_progress, advance=1)
                #for i in range(n_rows):
                # [3] Calculate the VFRAME (m/s) for each VEGAS data row
                vframe_i = calc_vframe(lock, gbt, session.DMJD, session.RA, session.DEC, session.GO_VELDEF)
                lo1a_dict[3]['data']['VFRAME'] = vframe_i
                # [4] Calculate the new RVSYS values
                rvsys_i = calc_rvsys(lock, session.LO1_VFRAME, vframe_i, session.LO1_RVSYS)
                lo1a_dict[3]['data']['RVSYS'] = rvsys_i
                # [5] Calculate the new LO1 frequencies
                lo1freq_i = calc_lo1freq(lock, session.LO1_RESTFRQ, rvsys_i, session.LO1_LOMULT, session.LO1_LOOFFSET, session.LO1_IFFREQ, session.LO1_SIDEBAND)
                lo1a_dict[3]['data']['LO1FREQ'] = lo1freq_i
                # [6] Calculate a frequency offset (F_OFFSET) from VFRAME
                f_sky_i, f_offset_i = calc_f_offset(lock, session.LO1_LO1FREQ, lo1freq_i, session.LO1_LOMULT, session.LO1_LOOFFSET, session.LO1_IFFREQ, session.LO1_SIDEBAND)
                # [7] Calculate a channel shift from the frequency offset
                channel_shift_i = calc_channel_offset(lock, f_offset_i, session.VEGAS_CDELT1, session.LO1_SIDEBAND)
                # [8] Shift the rows of data in the VEGAS FITS files by their channel offsets
                for csi in range(len(channel_shift_i)):
                    vegas_dict[6]['data']['data'][csi] = shift(vegas_dict[6]['data']['data'][csi], [0,0,channel_shift_i[csi]], cval=np.NaN)
                # [9] Write the DMJD, RA, and DEC values in the row of the LO1TBL
                lo1a_dict[3]['data']['DMJD'] = session.DMJD
                lo1a_dict[3]['data']['RA'] = session.RA
                lo1a_dict[3]['data']['DEC'] = session.DEC
                
                progress.update(task_id, advance=1)
                progress.update(overall_progress, advance=1)
                    
                # [10] Shift the other VEGAS values
                n_cri = np.shape(session.VEGAS_CRVAL1)[0]
                n_sp = np.shape(session.VEGAS_SPURCHAN)[0]
                for cri in range(n_cri):
                    vegas_dict[4]['data']['CRVAL1'][cri] = float(vegas_dict[4]['data']['CRVAL1'][cri] + f_offset_i[cri])

                for sp in range(n_sp):
                    vegas_dict[1]['data']['SPURCHAN'][sp] = int((sp*session.VEGAS_ADCSAMPF/64-vegas_dict[4]['data']['CRVAL1'][sp])/session.VEGAS_CDELT1[sp]+session.VEGAS_CRPIX1)
                    vegas_dict[1]['data']['SPURFREQ'][sp] = float(vegas_dict[1]['data']['SPURFREQ'][sp]+f_offset_i[sp])
                
                for vk in list(vegas_dict[0]['hdr']):
                    if vk.startswith('SUB'):
                        vegas_dict[0]['hdr'][vk] = float(vegas_dict[0]['hdr'][vk]+f_offset_i[0])
                progress.update(task_id, advance=1)
                progress.update(overall_progress, advance=1)

                # Writing the new VEGAS files
                write_new_VEGAS(v_path_og, v_path_save, vegas_dict)
                progress.update(task_id, advance=1)
                progress.update(overall_progress, advance=1)

                # Writing the new LO1A files
                write_new_LO1A(l_path_og, l_path_save, lo1a_dict)
                progress.update(task_id, advance=1)
                progress.update(overall_progress, advance=1)