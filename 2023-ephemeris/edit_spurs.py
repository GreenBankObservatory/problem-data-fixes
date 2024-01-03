from glob import glob
from astropy.io import fits
import os, multiprocessing
import numpy as np
from tqdm import tqdm
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from rich.progress import (
    BarColumn,
    TimeElapsedColumn,
    Progress,
    TaskID
)

def get_threads():
    """ Get the number of CPU threads """
    max_threads = multiprocessing.cpu_count()
    thread_count = max_threads
    return thread_count

def make_progress_bar():
    """ Make a progress bar to keep track of the results """
    progress_bar = Progress(
        "[progress.description]{task.description}",
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        TimeElapsedColumn(),
        refresh_per_second=1,
        transient=True,
    )
    return progress_bar

def main_fix(progress, task_id, session):
    progress.start_task(task_id)
    vs = glob(os.path.join(session, "VEGAS/*.fits"))
    for v in vs:
        v_og = v.replace("modified", "original")
        vdata = fits.open(v, cache=False, lazy_load_hdus=True)
        vdata_og = fits.open(v_og)
        spurs = vdata[1].data['SPURCHAN']
        n_nan = np.count_nonzero(np.isnan(vdata[6].data['DATA'][0][0][0]))
        spur_shift = int(n_nan - 1)
        for si, s in enumerate(spurs):
            vdata[1].data['SPURCHAN'][si] = vdata_og[1].data['SPURCHAN'][si] + spur_shift
        for i in range(8):
            try:
                subfreq = vdata[0].header[f'SUB{i}FREQ']
                if type(subfreq) == float:
                    subfreq = str(int(subfreq)) + "."
                vdata[0].header[f'SUB{i}FREQ'] = subfreq
            except:
                pass
        vdata.writeto(v, overwrite=True, output_verify='ignore')
        vdata.close()
        vdata_og.close()
        del vdata
        del vdata_og
        progress.update(task_id, advance=1)

if __name__ == "__main__":

    thread_count = get_threads()
    progress_bar = make_progress_bar()

    fpaths = glob("/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/*")

    with progress_bar:
        overall_progress = progress_bar.add_task(f"[green]Spur Fix Progress:")
        sessions_to_fix = sorted(fpaths)
        sessions_to_fix.reverse()

        if len(sessions_to_fix) > 0:
            task_ids = []
            for session in sessions_to_fix:
                n_vegas = len(glob(os.path.join(session, "VEGAS/*.fits")))
                task_id = progress_bar.add_task(f"[cyan]{session}", start=False, transient=True, total=n_vegas)
                task_ids.append(task_id)
            done_count = 0
            done_total = len(sessions_to_fix)
            progress_bar.update(overall_progress)

            lock = threading.Lock()
            with ThreadPoolExecutor(thread_count) as executor:
                for out in as_completed([executor.submit(main_fix, progress_bar,  task_id, session) 
                                        for session, task_id in zip(sessions_to_fix, task_ids)]):
                    result = out.result()
                    if result is not None:
                        print(result)
                        print("SESSION WITH ERROR: ", session)
                    done_count += 1
                    progress_bar.update(overall_progress, completed=done_count, total=done_total)