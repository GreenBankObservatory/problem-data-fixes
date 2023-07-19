import os, glob
import numpy as np
from ephemeris_fix import fix_lo1a_vegas

def welcome_message():
    n_char = os.get_terminal_size().columns
    message = "EPHEMERIS FIX"
    if n_char >= 25:
        full_line = "#"*n_char
        space_line = "#"*3 + " "*(n_char-6) + "#"*3
        lpad = int(np.floor((n_char-19)/2))
        rpad = int(np.ceil((n_char-19)/2))
        text_line = "#"*3 + " "*lpad + message + " "*rpad + "#"*3 
        print(full_line)
        print(full_line)
        print(space_line)
        print(text_line)
        print(space_line)
        print(full_line)
        print(full_line)
    else:
        print(message)

def check_host(fpath="."):
    return os.path.exists(fpath)

def check_project(project_code):
    fpath = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/original/{project_code}*"
    sessions = glob.glob(fpath)
    return os.path.exists(fpath), sessions

welcome_message()

data_path = "/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/original/"
save_data_path = "/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/"
avail_sessions = [g.split('/')[-1] for g in glob.glob(f"{data_path}*")]
print("Please select from the following sessions: ")
print(avail_sessions)
selected_session = input("Session: ")
if check_host(data_path + selected_session):
    print(f"Great! I'll get started on {selected_session}.")
    fix_lo1a_vegas(data_path, selected_session)
    print(f"All done! You can find the data at {save_data_path}{selected_session}.")
else:
    print("Sorry, I couldn't find that session. Please try again.")