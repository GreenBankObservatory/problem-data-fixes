import pandas as pd
import numpy as np
from util.core import *

def read_data():
    true_values = pd.read_csv("./tests/local-data/true_values.csv")
    return true_values

def get_vframe_error(veldef):
    gbt = create_obj_gbt()
    true_vals = read_data()
    indx = true_vals['VELDEF'] == veldef
    vframes = calc_vframe(gbt, true_vals['DMJD'][indx], true_vals['RA'][indx], true_vals['DEC'][indx], veldef)
    vframe_diffs = np.abs(np.array(true_vals['VFRAME'][indx]) - np.array(vframes))
    return np.max(vframe_diffs)

def test_vframe_barycentric():
    bary_err = get_vframe_error('VOPT-BAR')
    assert bary_err <= 1.0

def test_vframe_heliocentric():
    helio_err = get_vframe_error('VOPT-HEL')
    assert helio_err <= 1.0

def test_vframe_LSRK():
    lsrk_err = get_vframe_error('VRAD-LSR')
    assert lsrk_err <= 1.0