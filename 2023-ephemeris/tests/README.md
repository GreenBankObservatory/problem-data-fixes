# 2023-ephemeris Tests
These are tests that run for the 2023-ephemeris code. They are as follows:

## test_fits_writer.py
* Criteria: The old FITS file must remain unaffected, and the new file must have a new FITS keyword. 

## test_lo1freq_values.py
* Criteria: LO1FREQ values must be within 150 Hz of the true values when calculated for test data

## test_vframe_values.py
* Criteria: VFRAME values must be within 1 m/s of the true values when calculated for test data