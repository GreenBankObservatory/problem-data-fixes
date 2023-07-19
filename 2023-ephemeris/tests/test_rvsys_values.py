import pandas as pd
import numpy as np
from core import *

def read_data():
    true_values = pd.read_csv("./tests/data/true_values_new.csv")
    return true_values

def test_rvsys():
    data = read_data()
    rvsys = calc_rvsys(data['VEL'], data['VFRAME'])
    rvsys_err = np.abs(np.subtract(data['RVSYS'], rvsys))
    assert np.max(rvsys_err) <= 5000