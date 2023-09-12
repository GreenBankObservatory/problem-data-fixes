import os, glob
import numpy as np
from ephemeris_fix import *
from concurrent.futures import ThreadPoolExecutor
from rich.progress import (
    BarColumn,
    DownloadColumn,
    TextColumn,
    TransferSpeedColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
    Progress,
    TaskID,
)

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
all_sessions = [g.split('/')[-1] for g in glob.glob(f"{data_path}*")]
avail_sessions = [a for a in all_sessions if a .startswith("AGBT")]
#print(avail_sessions)
done_sessions = ["AGBT21B_316_21", "AGBT23A_387_01", "AGBT23A_387_02", "AGBT23A_387_03", "AGBT23A_387_04"]
#print("Please select from the following sessions: ")
#print(avail_sessions)
#selected_session = "AGBT21B_316_21"#input("Session: ")
# selected_session = input("Session: ")# "AGBT23A_387_01"
##s_tot = 10
##sn = int(input("Skip number: "))
#n_avail = len(avail_sessions)
#n_done = 0
#for selected_session in avail_sessions[sn::s_tot]:
#    if selected_session not in done_sessions:
#for selected_session in done_sessions:
#    if check_host(data_path + selected_session):
#        print(f"Great! I'll get started on {selected_session}.")
#        fix_lo1a_vegas(data_path, selected_session)
#        #fix_sdfits(selected_session)
#        print(f"All done! You can find the data at {save_data_path}{selected_session}.")
#    else:
#        print("Sorry, I couldn't find that session. Please try again.")

n_processes = 12

progress_bar = Progress(
    "[progress.description]{task.description}",
    BarColumn(),
    "[progress.percentage]{task.percentage:>3.0f}%",
    TimeRemainingColumn(),
    TimeElapsedColumn(),
    refresh_per_second=1,
)

with progress_bar:
    # Create a task for the progress bar
    overall_progress = progress_bar.add_task(f"[green]Ephemeris:")

    # Define the processes
    futures = []
    executor = ThreadPoolExecutor(n_processes)
    for avail_session in avail_sessions:
        futures.append(executor.submit(fix_lo1a_vegas, data_path, avail_session))

    # Update the progress bar
    while (n_finished := sum([future.done() for future in futures])) < len(futures):
        progress_bar.update(overall_progress, completed=n_finished, total=len(futures))