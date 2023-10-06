from glob import glob
import shutil, os
from tqdm import tqdm

fix_path = "/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/*"
session_codes = glob(fix_path)
ignore_dir = ['LO1A', 'VEGAS', 'VEGAS_CODD', 'ScanLog.fits']
for sc_path in session_codes:
    sc = sc_path.split('/')[-1]
    fix_dirs = glob(sc_path + '/*')
    old_dirs = glob(f"/home/gbtdata/{sc}/*")

    print(sc)
    print('--------------------------')
    for od_path in old_dirs:
        od = od_path.split('/')[-1]
        if od not in ignore_dir:
            od_fits = glob(od_path + '/*.fits*')
            print(od, len(od_fits))
            os.mkdir(sc_path + '/' + od)
            for odf in tqdm(od_fits):
                #print(f"Copying {odf} to {sc_path}/{od}/")
                shutil.copy(odf, sc_path + '/' + od + '/')
        #if od == 'ScanLog.fits':
        #    print(f"Moving {od_path} to {sc_path}/")
        #    shutil.copy(od_path, sc_path + '/')
    print('\n')