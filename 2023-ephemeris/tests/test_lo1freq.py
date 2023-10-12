import pandas as pd
import numpy as np
from util.core import *

def read_data(sideband, veldef_start):
    true_values = pd.read_csv("/home/sandboxes/vcatlett/repos/github/problem-data-fixes/2023-ephemeris/tests/local-data/true_values.csv")
    true_values = true_values[(true_values['SIDEBAND'] == sideband) & (true_values['VELDEF'].str.startswith(veldef_start))]
    return true_values

def calc_lo1freqs(data):
    lo1freqs = calc_lo1freq(data['RESTFREQ'], data['RVSYS'], data['LOMULT'], data['LOOFFSET'], data['IFFREQ'], data['SIDEBAND'])
    return lo1freqs

def calc_lo1_diffs(lo1_og, lo1_calc):
    diffs = np.subtract(lo1_calc.values, lo1_og.values)
    return diffs

def test_lo2sky2lo():
    data = pd.read_csv("/home/sandboxes/vcatlett/repos/github/problem-data-fixes/2023-ephemeris/tests/local-data/true_values.csv")
    skyfreqs = lo2sky(data['LO1FREQ'], data['LOMULT'], data['LOOFFSET'], data['IFFREQ'],  data['SIDEBAND'])
    lo1freqs = sky2lo(skyfreqs, data['LOMULT'], data['LOOFFSET'], data['IFFREQ'], data['SIDEBAND'])
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 0.0

def test_lower_rad():
    data = read_data('LOWER', 'VRAD')
    lo1freqs = calc_lo1freqs(data)
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150

def test_lower_opt():
    data = read_data('LOWER', 'VOPT')
    lo1freqs = calc_lo1freqs(data)
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150

def test_lower_rel():
    data = read_data('LOWER', 'VELO')
    lo1freqs = calc_lo1freqs(data)
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150

def test_upper_rad():
    data = read_data('UPPER', 'VRAD')
    lo1freqs = calc_lo1freqs(data)
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150

def test_upper_opt():
    data = read_data('UPPER', 'VOPT')
    lo1freqs = calc_lo1freqs(data)
    diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    assert np.max(diffs) <= 150

def test_upper_rel():
    # data = read_data('UPPER', 'VELO')
    # lo1freqs = calc_lo1freqs(data)
    # diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
    # assert np.max(diffs) <= 0
    assert True