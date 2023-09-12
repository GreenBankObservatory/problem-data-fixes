import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from core import *

def read_data(sideband, veldef_start):
    true_values = pd.read_csv("./tests/local-data/true_values.csv")
    true_values = true_values[(true_values['SIDEBAND'] == sideband) & (true_values['VELDEF'].str.startswith(veldef_start))]
    return true_values

def calc_lo1freqs(data, sideband, vel_formula):
    lo1freqs = sky2lo(data['RESTFREQ'], data['LOMULT'], data['IFFREQ'], data['VFRAME'], sideband, vel_formula, data['VEL'])
    return lo1freqs

def calc_lo1freqs_prop(data, sideband, vel_formula):
    lo1freqs = sky2lo(data['RESTFREQ'], data['LOMULT'], data['IFFREQ'], data['VFRAME_CALC'], sideband, vel_formula, data['VEL'])
    return lo1freqs

def calc_lo1_diffs(lo1_og, lo1_calc):
    diffs = np.abs(np.subtract(lo1_calc, lo1_og))
    return diffs

def read_vframe_data():
    true_values = pd.read_csv("./tests/local-data/true_values.csv")
    return true_values

def get_vframe_error(veldef):
    true_vals = read_vframe_data()
    indx = true_vals['VELDEF'] == veldef
    vframes = calc_vframe(true_vals['DMJD'][indx], true_vals['RA'][indx], true_vals['DEC'][indx], veldef)
    vframe_diffs = np.abs(np.array(true_vals['VFRAME'][indx]) - np.array(vframes))
    return vframes, vframe_diffs

true_values = pd.read_csv("./tests/local-data/true_values.csv")

#rvsys = calc_rvsys(true_values['VEL'], true_values['VFRAME'])

print("VFRAME...")
all_vframes = []
all_vframe_err = []
for vd in ['VOPT-BAR', 'VOPT-HEL', 'VRAD-LSR']:
    vframes, vframe_err = get_vframe_error(vd)
    all_vframe_err = all_vframe_err + list(vframe_err)
    all_vframes = all_vframes + list(vframes)
#true_values['VFRAME_CALC'] = all_vframes

print("LO1FREQ...")
all_lo1_diffs = []
all_lo1_prop_diffs = []
for sb in ["LOWER", "UPPER"]:
    for vd in ['VOPT', 'VRAD']:
        data = true_values[(true_values['SIDEBAND'] == sb) & (true_values['VELDEF'].str.startswith(vd))]
        lo1freqs = calc_lo1freqs(data, sb, vd[1:].lower())
        diffs = calc_lo1_diffs(data['LO1FREQ'], lo1freqs)
        all_lo1_diffs = all_lo1_diffs + list(diffs)
        
        #lo1freqs_prop = calc_lo1freqs_prop(data, sb, vd[1:].lower())
        #diffs_prop = calc_lo1_diffs(data['LO1FREQ'], lo1freqs_prop)
        #all_lo1_prop_diffs = all_lo1_prop_diffs + list(diffs_prop)

print("RVSYS...")
rvsys= calc_rvsys(true_values['VEL'], true_values['VFRAME'], true_values['VELDEF'])
rvsys_err = np.abs(np.subtract(true_values['RVSYS'], rvsys))
#rvsys_prop = calc_rvsys(true_values['VEL'], true_values['VFRAME_CALC'], true_values['VELDEF'])
#rvsys_prop_err = np.abs(np.subtract(true_values['RVSYS'], rvsys_prop))

fig, ax = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(15, 5))

cax = ax[0]
cax.set_title("VFRAME")
cax.hist(all_vframe_err, bins=10)
cax.set_xlabel("| VFRAME Error | (m/s)")
cax.set_ylabel("Counts")

cax = ax[1]
cax.set_title("LO1FREQ")
cax.hist(all_lo1_diffs, bins=10)
cax.set_xlabel("| LO1FREQ Error | (Hz)")
#cax.set_ylabel("Counts")

cax = ax[2]
cax.set_title("RVSYS")
cax.hist(rvsys_err, bins=10)
cax.set_xlabel("| RVSYS Error | (m/s)")
#cax.set_ylabel("Counts")
'''
cax = ax[1][1]
cax.set_title("LO1FREQ")
cax.hist(all_lo1_prop_diffs, bins=10)
cax.set_xlabel("| LO1FREQ Error | (Hz)")
cax.set_ylabel("Counts")

cax = ax[1][2]
cax.set_title("RVSYS")
cax.hist(rvsys_prop_err, bins=10)
cax.set_xlabel("| RVSYS Error | (m/s)")
#cax.set_ylabel("Counts")
'''
plt.show()
plt.savefig("/home/sandboxes/vcatlett/repos/github/problem-data-fixes/2023-ephemeris/EPHEMERIS_ERRORS.png")
plt.close()
#fig, ax = plt.subplots(nrows=2, ncols=3)
#xs = ['DMJD', 'VEL', 'VFRAME']

'''
caxi = 0
cax = ax[0][caxi]
cax.scatter(true_values[xs[caxi]], true_values['RVSYS'], color='k')
cax.scatter(true_values[xs[caxi]], rvsys, color='r')

cax = ax[1][caxi]
cax.scatter(true_values[xs[caxi]], rvsys_err, color='r')
cax.set_xlabel(xs[caxi])

caxi = 1
cax = ax[0][caxi]
cax.scatter(true_values[xs[caxi]], true_values['RVSYS'], color='k')
cax.scatter(true_values[xs[caxi]], rvsys, color='r')

cax = ax[1][caxi]
cax.scatter(true_values[xs[caxi]], rvsys_err, color='r')
cax.set_xlabel(xs[caxi])

caxi = 2
cax = ax[0][caxi]
cax.scatter(true_values[xs[caxi]], true_values['RVSYS'], color='k')
cax.scatter(true_values[xs[caxi]], rvsys, color='r')

cax = ax[1][caxi]
cax.scatter(true_values[xs[caxi]], rvsys_err, color='r')
cax.set_xlabel(xs[caxi])

plt.show()
'''