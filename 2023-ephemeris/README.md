# 2023-ephemeris
This fix is for data which had improper Doppler tracking after the JPL ephemeris file silently expired. This encompasses VEGAS spectral line data taken in non-topocentric reference frames between MJD 60016 at 00:00 UTC and MJD 60090 at 21:01 UTC. 

## Set up

(1) ```$ git clone git@github.com:GreenBankObservatory/problem-data-fixes.git```

(2) ```$ cd problem-data-fixes/2023-ephemeris```

(3) ```$ /users/gbosdd/python/bin/python3.11 -m venv ephem-env-3.11```

(4) ```$ source ephem-env-3.11/bin/activate```

(5) ```(ephem-env-3.11) $ pip install -r requirements.txt```

(6) ```(ephem-env-3.11) $ python main.py```