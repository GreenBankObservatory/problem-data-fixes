import pandas as pd
import numpy as np
from util.core import *

def read_data():
    true_values = pd.read_csv("/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/testing/pytest_data/true_values.csv")
    return true_values

def get_rvsys_error_frame(veldef):
    lock = None
    data = read_data()
    indx = data['VELDEF'].str.endswith(veldef)
    rvsys = calc_rvsys(lock, data['VFRAME'][indx], data['VFRAME'][indx], data['RVSYS'][indx])
    rvsys_err = np.abs(np.subtract(data['RVSYS'][indx].values, rvsys))
    return np.max(rvsys_err)

def get_rvsys_error_fdef(fdef):
    lock = None
    data = read_data()
    indx = data['VELDEF'].str.startswith(fdef)
    rvsys = calc_rvsys(lock, data['VFRAME'][indx], data['VFRAME'][indx], data['RVSYS'][indx])
    rvsys_err = np.abs(np.subtract(data['RVSYS'][indx].values, rvsys))
    return np.max(rvsys_err)
    
def test_rvsys_barycentric():
    """ Check that max barycentric RVSYS error is <= 0.05 m/s """
    rvsys_err = get_rvsys_error_frame('BAR')
    assert np.max(rvsys_err) <= 1E-6
    
def test_rvsys_heliocentric():
    """ Check that max heliocentric RVSYS error is <= 0.05 m/s """
    rvsys_err = get_rvsys_error_frame('HEL')
    assert np.max(rvsys_err) <= 1E-6

def test_rvsys_LSR():
    """ Check that max LSR RVSYS error is <= 0.05 m/s """
    rvsys_err = get_rvsys_error_frame('LSR')
    assert np.max(rvsys_err) <= 1E-6

def test_rvsys_topocentric():
    """ Check that max topocentric RVSYS error is = 0.0 """
    rvsys_err = get_rvsys_error_frame('TOP')
    assert np.max(rvsys_err) <= 1E-6

def test_rvsys_radio():
    """ Check that max radio RVSYS error is <= 0.05 m/s """
    rvsys_err = get_rvsys_error_fdef('VRAD')
    assert np.max(rvsys_err) <= 1E-6

def test_rvsys_opt():
    """ Check that max optical RVSYS error is <= 0.05 m/s """
    rvsys_err = get_rvsys_error_fdef('VOPT')
    assert np.max(rvsys_err) <= 1E-6

def test_rvsys_rel():
    """ Check that max relativistic RVSYS error is <= 0.05 m/s """
    rvsys_err = get_rvsys_error_fdef('VELO')
    assert np.max(rvsys_err) <= 1E-6