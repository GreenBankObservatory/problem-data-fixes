import pandas as pd
import numpy as np
from util.core import *

def read_data():
    """ Open the test data """
    true_values = pd.read_csv("/home/sandboxes/vcatlett/repos/github/problem-data-fixes/2023-ephemeris/tests/local-data/true_values.csv")
    return true_values

def get_vframe_error(veldef):
    """ Run the test for a given VELDEF """
    gbt = create_obj_gbt()
    true_vals = read_data()
    indx = true_vals['VELDEF'].str.endswith(veldef)
    vframes = calc_vframe(gbt, true_vals['DMJD'][indx], true_vals['RA'][indx], true_vals['DEC'][indx], true_vals['VELDEF'][indx])
    vframe_diffs = np.abs(np.subtract(true_vals['VFRAME'][indx], vframes))
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