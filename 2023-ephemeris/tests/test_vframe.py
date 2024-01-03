import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from util.core import *

def read_data():
    """ Open the test data """
    true_values = pd.read_csv("/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/pytest_data/true_values.csv")
    return true_values

def get_vframe_error(veldef):
    """ Run the test for a given VELDEF """
    lock = None
    gbt = create_obj_gbt(lock)
    true_vals = read_data()
    indx = true_vals['VELDEF'].str.endswith(veldef)

    n_max = 10000 # To run faster and prevent memory errors
    true_vals_small = true_vals[indx]
    if len(true_vals_small) > n_max:
        true_vals_small = true_vals_small[0:n_max]
    true_vals_small = true_vals_small.reset_index()

    vframes = calc_vframe(lock, gbt, true_vals_small['DMJD'], true_vals_small['RA'], true_vals_small['DEC'], true_vals_small['VELDEF'])
    vframe_diffs = np.abs(np.subtract(true_vals_small['VFRAME'], vframes))
    plt.hist(vframe_diffs)
    plt.title(f"Errors for VFRAME = {veldef} (n = {len(true_vals_small)})")
    plt.xlabel("VFRAME Error (m/s)")
    plt.savefig(f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/pytest_data/plots/ERRORS-VFRAME-{veldef}.png")
    return np.max(vframe_diffs)

def test_vframe_barycentric():
    """ Check that max barycentric VFRAME error is <= 0.05 m/s """
    bary_err = get_vframe_error('BAR')
    assert bary_err <= 0.05

def test_vframe_heliocentric():
    """ Check that max heliocentric VFRAME error is <= 0.05 m/s """
    helio_err = get_vframe_error('HEL')
    assert helio_err <= 0.05

def test_vframe_LSR():
    """ Check that max LSR VFRAME error is <= 0.05 m/s """
    lsrk_err = get_vframe_error('LSR')
    assert lsrk_err <= 0.05

def test_vframe_topocentric():
    """ Check that max topocentric VFRAME error is = 0.0 """
    lsrk_err = get_vframe_error('TOP')
    assert lsrk_err <= 0.0