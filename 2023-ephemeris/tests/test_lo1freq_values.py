import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from util.core import *

def read_data(sideband, veldef_start):
    true_values = pd.read_csv("./tests/local-data/true_values.csv")
    true_values = true_values[(true_values['SIDEBAND'] == sideband) & (true_values['VELDEF'].str.startswith(veldef_start))]
    return true_values

def calc_lo1freqs(data, sideband, vel_formula):
    f_new, f_offset = calc_f_offset(data['RESTFREQ'], data['VFRAME'], data['LO1FREQ'], data['LOMULT'], data['IFFREQ'], sideband, formula=vel_formula)
    lo1freqs = sky2lo(f_new, data['LOMULT'], data['LOOFFSET'], data['IFFREQ'], sideband)
    return lo1freqs

def calc_lo1_diffs(lo1_og, lo1_calc):
    diffs = np.abs(np.subtract(lo1_calc, lo1_og))
    return diffs

def test_lower_rad():
    data = read_data('LOWER', 'VRAD')
    lo1freqs = calc_lo1freqs(data, 'LOWER', 'rad')
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150

def test_lower_opt():
    data = read_data('LOWER', 'VOPT')
    lo1freqs = calc_lo1freqs(data, 'LOWER', 'opt')
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150

def test_upper_rad():
    data = read_data('UPPER', 'VRAD')
    lo1freqs = calc_lo1freqs(data, 'UPPER', 'rad')
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150

def test_upper_opt():
    data = read_data('UPPER', 'VOPT')
    lo1freqs = calc_lo1freqs(data, 'UPPER', 'opt')
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150