import glob
import pandas as pd
from astropy.io import fits
from datetime import datetime
from tqdm import tqdm
import numpy as np

def get_uid(session):
    session_split = session.split("_")
    session_uid = ''
    if session_split[-1].isnumeric():
        n_concat = len(session_split)-1
        session_num = session_split[-1]
    else:
        n_concat = len(session_split)
        session_num = "0"
    for i in range(n_concat):
        if i == 0:
            session_uid = session_split[i]
        else:
            session_uid = f"{session_uid}_{session_split[i]}"
    return session_uid, session_num

def load_data():
    save_dict = {}

    data_path = "/home/gbtdata/"
    all_sessions = [g.split('/')[-1] for g in glob.glob(f"{data_path}*")]
    date_format = '%Y-%m-%dT%H:%M:%S'
    t_start = datetime(2023, 3, 13, 0, 0, 0)
    t_end = datetime(2023, 5, 26, 21, 0, 1)
    
    unique_sessions = []
    for a in tqdm(all_sessions):        
        session_uid = get_uid(a)[0]
        if session_uid not in unique_sessions:
            unique_sessions.append(session_uid)

    for a in tqdm(all_sessions):
        go_files = glob.glob(f"/home/gbtdata/{a}/GO/*.fits")
        project_id_both = get_uid(a)
        project_id = project_id_both[0]
        session_id = int(project_id_both[1])
        a_instrument = []
        a_frame = []
        a_observer = []
        if a.startswith("AGBT1B_316"):
            print(session_id)
        for g in go_files:
            g_data = fits.open(g)
            dt_o_og = g_data[0].header['DATE-OBS']
            dt_o = datetime.strptime(dt_o_og, date_format)
            i_o_all = g_data[0].header['HISTORY']
            for ioi in i_o_all:
                if ioi.startswith("backend"):
                    i_o = ioi.split("=")[1].replace("\'", '').replace(" ", '')  
            #breakpoint()
            f_o = g_data[0].header['VELDEF']
            a_o = g_data[0].header['OBSERVER']

            if (t_start <= dt_o) and (dt_o <= t_end):
                if i_o not in a_instrument:
                    a_instrument.append(i_o)
                if f_o not in a_frame:
                    a_frame.append(f_o)
                if a_o not in a_observer:
                    a_observer.append(a_o)
            
            g_data.close()

        if project_id not in save_dict.keys():
            if a.startswith("AGBT1B_316"):
                print("Starting to save...")
            save_dict[project_id] = {
                'session': [session_id],
                'instrument':a_instrument,
                'frame':a_frame,
                'observer':a_observer
                }
        else:
            save_dict[project_id]['session'].append(session_id)
            for f in a_frame:
                if f not in save_dict[project_id]['frame']:
                    save_dict[project_id]['frame'].append(f)
            for o in a_observer:
                if o not in save_dict[project_id]['observer']:
                    save_dict[project_id]['observer'].append(o)
            for i in a_instrument:
                if i not in save_dict[project_id]['instrument']:
                    save_dict[project_id]['instrument'].append(i)

    save_df_dict = {
        'project':[],
        'session':[],
        'instrument':[],
        'frame':[],
        'observer':[]
    }
    for p in save_dict.keys():
        pdict = save_dict[p]
        save_df_dict['project'].append(p)
        save_df_dict['session'].append(pdict['session'])
        save_df_dict['instrument'].append(pdict['instrument'])
        save_df_dict['frame'].append(pdict['frame'])
        save_df_dict['observer'].append(pdict['observer'])
    df = pd.DataFrame(save_df_dict)
    df.to_csv('static/All_Sessions.csv', index=False)

def make_range_str(session_list):
    session_list.sort()
    max_i = len(session_list)-1
    list_str = ''
    for i in range(len(session_list)):
        si = session_list[i]
        
        if i == 0:
            list_str = str(si)
            cont_bool = True
        elif i == max_i:
            list_str = f"{list_str}-{si}"
        else:
            if cont_bool:
                if si + 1 == session_list[i+1]:
                    cont_bool=True
                else:
                    list_str = f"{list_str}-{si}"
                    cont_bool=False
            else:
                list_str = f"{list_str}, {si}"
                if si + 1 == session_list[i+1]:
                    cont_bool=True
                
    return list_str

def create_table_string(highlight=(False, False, False),table_class='wikitable', style=''):
    """
    Takes a list and returns a wikitable.
    @param data: The list that is converted to a wikitable.
    @type data: List (Nested)
    @param highlight: Tuple of rows and columns that should be highlighted.
                      (first row, last row, left column, right column)
    @type highlight: Tuple
    @param table_class: A string containing the class description.
                        See wikitable help.
    @type table_class: String
    @param style: A string containing the style description.
                  See wikitable help.
    @type style: String
    """
    data_df = pd.read_csv("static/All_Sessions.csv")
    data_df = data_df.sort_values(by='project')
    data_df.reset_index()
    data_project = list(data_df["project"])
    data_session = list(data_df["session"])
    data_instrument = list(data_df["instrument"])
    data_frame = list(data_df["frame"])
    data_observer = list(data_df["observer"])
    data = [
        [
            [ data_project[i], data_session[i], data_instrument[i], data_frame[i], data_observer[i] ]
            for i in range(len(data_session))
        ]
    ]
    data = data[0]
    my_name_str = "[[Main.VictoriaCatlett][Victoria Catlett]]"
    no_fix_str = ""
    with open('static/wikitable_all.txt', 'w') as the_file:
        header_str = "| Project | Session(s) | Affected? | Backend(s) | VELDEF | Observer | Project Friend | Approved? | Approved By | Comments |\n"
        the_file.write(header_str)
        for key, row in tqdm(enumerate(data)):
            # projects,session,frame,observer
            proj = row[0]
            sesh = row[1].replace("[", '').replace("]", '').replace("\'", '').replace("\"", '')
            inst = row[2].replace("[", '').replace("]", '').replace("\'", '').replace("\"", '')
            veldefs = row[3].replace("[", '').replace("]", '').replace("\'", '').replace("\"", '')
            obs = row[4].replace("[", '').replace("]", '').replace("\'", '').replace("\"", '')
            sesh_split = sesh.split(',')
            sesh_split = [int(s) for s in sesh_split]
            if len(sesh_split)>1:
                sesh_split.sort()
                sesh_str = make_range_str(sesh_split)
            elif int(sesh[0]) == 0:
                sesh_str = "N/A"
            else:
                sesh_str=str(sesh[0])
                
            if proj.startswith("TRFI"):
                veldefs = "VRAD-TOP"
            need_fix = True
            if len(veldefs) == 8: # Only 1 veldef
                if veldefs.endswith("TOP"):
                    need_fix = False
                    no_fix_str = "No fix needed (topocentric)"
            all_backend = inst.replace(" ", '').split(',')
            #print(len(all_backend), len(all_backend[0]))
            if len(all_backend[0]) > 0:
                if "VEGAS" not in all_backend:
                    need_fix = False
                    no_fix_str = "No fix needed (not VEGAS)"

            if need_fix:
                row_str = f'| {proj} | {sesh_str} | Yes | {inst} | {veldefs} | {obs} | - | %X% | - | - |\n'
            else:
                row_str = f'| {proj} | {sesh_str} | No | {inst} | {veldefs} | {obs} | - | %Y% | {my_name_str} | {no_fix_str} |\n'
            the_file.write(row_str)

if __name__ == "__main__":
    #load_data()
    create_table_string()