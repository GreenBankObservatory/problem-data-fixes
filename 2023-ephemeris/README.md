# 2023-ephemeris
This fix is for data which had improper Doppler tracking after the JPL ephemeris file silently expired. This encompasses VEGAS spectral line data taken in non-topocentric reference frames between MJD 60016 at 00:00 UTC and MJD 60090 at 21:01 UTC. 

## Set up

(1) ``$ git clone git@github.com:GreenBankObservatory/problem-data-fixes.git``

(2) ``$ cd problem-data-fixes/2023-ephemeris``

(3) ``$ /users/gbosdd/python/bin/python3.11 -m venv ephem-env-3.11``

(4) ``$ source ephem-env-3.11/bin/activate``

(5) ``(ephem-env-3.11) $ pip install -r requirements.txt``

(6) ``(ephem-env-3.11) $ python main.py``

(7) Follow the prompts. Use glob-like selection for choosing sessions (ex. ``AGBT*``)

To check on the overall memory usage of the script, log into the same machine in another terminal. Then run ``ps aux --sort=-%mem | grep {USERNAME} | head``. If using more than 1 CPU thread, you will likely see a `%CPU` greater than 100. That's fine. Every 100% translates to full utilization of 1 core, so on a machine with 16 cores, the max utilization would be 1600%. I haven't seen this code get above 500% even with all 16 possible threads. 