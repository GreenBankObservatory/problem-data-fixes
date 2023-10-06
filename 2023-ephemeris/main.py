import os, glob, multiprocessing
from concurrent.futures import ThreadPoolExecutor
from getpass import getuser

from rich import print as rprint
from rich.prompt import Prompt
from rich.console import Console
from rich.table import Table
from rich.progress import (
    BarColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
    Progress,
    TaskID
)

from util.messages import Messages
from util.ephemeris_fix import *

def get_paths():
    """ Get paths to find the data and save the results """
    loadpath = get_loadpath()
    savepath = get_savepath()
    rprint(f"Loading data from {loadpath}")
    rprint(f"Saving results to {savepath}")
    return loadpath, savepath

def get_loadpath():
    """ Get path to find the data """
    loadpath_default = "/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/original"
    loadpath = Prompt.ask(f"Where can I find the files?", default = loadpath_default)
    loadpath_exists = check_path(loadpath)
    if loadpath_exists:
        return loadpath
    else:
        rprint(f"Sorry, I couldn't find that path. Please try again.")
        loadpath = get_loadpath()
        return loadpath

def get_savepath():
    """ Get path to save the results """
    USERNAME = getuser()
    if USERNAME == "vcatlett":
        savepath_default = "/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified"
    else:
        savepath_temp = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/{USERNAME}"
        if os.path.exists(savepath_temp):
            savepath_default = savepath_temp
        else:
            savepath_default = "."
    savepath = Prompt.ask("Where should I save the results?", default = savepath_default)
    savepath_exists = check_path(savepath)
    if savepath_exists:
        return savepath
    else:
        rprint(f"Sorry, I couldn't find that path. Please try again.")
        savepath = get_savepath()
        return savepath

def check_path(fpath):
    """ Check that the path exists """
    return os.path.exists(fpath)

def get_projects(loadpath):
    """ Find all projects on the loadpath """
    rprint("Searching for projects...")
    sessions = [ f.name for f in os.scandir(loadpath) if f.is_dir() ]
    projects = {}
    for s in sessions:
        s_split = s.split('_')
        if len(s_split) == 3:
            p = "_".join(s_split[0:2])
            sc = s_split[2]
        elif len(s_split) == 2:
            p = s_split[0]
            sc = s_split[1]
        if p not in projects.keys():
            projects[p] = [sc]
        else:
            projects[p].append(sc)
    get_project_table(projects)

def get_project_table(projects):
    """ Print a table of available projects """
    table = Table(title="Available Projects")
    table.add_column("Project", justify="right", style="cyan", no_wrap=True)
    table.add_column("Session(s)", justify="left", style="magenta")

    for p in sorted(projects.keys()):
        s_list = projects[p]
        if len(s_list) > 1:
            s_list = sorted(s_list)
            s_str = ", ".join(s_list)
        else:
            s_str = s_list[0]
        table.add_row(p, s_str)

    console = Console()
    console.print(table)

def select_sessions(loadpath):
    """ Get a list of sessions to fix """
    sessions_request = Prompt.ask(f"What sessions should I fix?", default = "*")
    session_paths = glob.glob(os.path.join(loadpath, sessions_request))
    sessions_to_fix = [os.path.basename(sp) for sp in session_paths]
    return sessions_to_fix

def get_threads():
    """ Get the number of CPU threads """
    max_threads = multiprocessing.cpu_count()
    thread_count = int(Prompt.ask("Number of CPU threads", default = str(max_threads)))
    if thread_count > max_threads:
        rprint(f"WARNING: {thread_count} is too many CPU threads. Defaulting to maximum value of {max_threads}.")
        thread_count = max_threads
    return thread_count

def make_progress_bar():
    """ Make a progress bar to keep track of the results """
    progress_bar = Progress(
        "[progress.description]{task.description}",
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        TimeElapsedColumn(),
        refresh_per_second=1
    )
    return progress_bar

def run_fix(loadpath, savepath, sessions_to_fix):
    thread_count = get_threads()
    progress_bar = make_progress_bar()
    print('\n')
    with progress_bar:
        overall_progress = progress_bar.add_task(f"[green]Ephemeris Fix Progress:")
        fix_processes = []
        executor = ThreadPoolExecutor(thread_count)
        for session in sorted(sessions_to_fix):
            task_id = progress_bar.add_task(f"[cyan]{session}", start=False, transient=True)
            future = executor.submit(main_fix, progress_bar, task_id, loadpath, savepath, session)
            fix_processes.append(future)
        while (n_finished := sum([fp.done() for fp in fix_processes])) < len(fix_processes):
            progress_bar.update(overall_progress, completed=n_finished, total=len(fix_processes))

if __name__ == "__main__":
    Messages.hello()
    loadpath, savepath = get_paths()
    get_projects(loadpath)
    sessions_to_fix = select_sessions(loadpath)
    run_fix(loadpath, savepath, sessions_to_fix)
    Messages.goodbye()