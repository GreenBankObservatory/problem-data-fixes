from astropy.io import fits
import numpy as np
from glob import glob
import os

session = "AGBT21B_316_24"
path_old = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/test_sessions/original/{session}"
path_new = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/test_sessions/modified/{session}"

vegas_fpaths = glob(f"{path_old}/VEGAS/*.fits*")
vegas_fnames = [os.path.basename(g) for g in vegas_fpaths]

comparisons = {
    'SPURCHAN': {
        'is_data': True,
        'hnum': 1,
    },
    'SPURFREQ': {
        'is_data': True,
        'hnum': 1,
    },
    'DATA': {
        'is_data': True,
        'hnum': 6,
    }
}

for v in vegas_fnames:
    print('\n', v)
    v_old = fits.open(os.path.join(path_old, 'VEGAS', v))
    v_new = fits.open(os.path.join(path_new, 'VEGAS', v))

    for k in comparisons.keys():
        hnum = comparisons[k]['hnum']
        if comparisons[k]['is_data']:
            vals_old = v_old[hnum].data[k]
            vals_new = v_new[hnum].data[k]
            if len(vals_old) <= 0:
                print(k, "No Old Data")
            elif len(vals_new) <= 0:
                print(k, "No New Data")
            else:
                vals_diff = np.max(np.abs(np.subtract(vals_new, vals_old)))
                if (vals_diff > 0):
                    print(k, vals_diff)

    v_old.close()
    v_new.close()