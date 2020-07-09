#!/usr/bin/env python
"""Utilities for generating alternate downtime schedules
"""

# imports
from os import path
from os import environ
import numpy as np
import pandas as pd
import sqlite3
import astropy
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import get_moon, get_sun, EarthLocation
from lsst.sims.utils import Site
from scipy.signal import argrelextrema


# constants
SITE = Site('LSST')
LOCATION = EarthLocation(lat=SITE.latitude, lon=SITE.longitude,height=SITE.height)
START_MJD = 59853
SURVEY_DURATION = 366*20
SYNODIC_MONTH_DAYS = 29.530589 # From Lang's *Astrophysical Data*, p. 57
DOWNTIME_DURATION = 14

REF_DOWNTIME_DB = path.join(environ['SIMS_DOWNTIMEMODEL_DIR'],
                            'data', 'scheduled_downtime.db')

def compute_full_moon_nights(location, start_mjd, duration=SURVEY_DURATION):
    # Use UTC midnight, shifted to the site using the longitude
    # Should be good to within ~15 minutes
    mjds = np.arange(start_mjd, start_mjd + duration) - LOCATION.lon.deg/360
    times = Time(mjds, format='mjd', location=LOCATION)
    
    sun_coords = get_sun(times)
    moon_coords = get_moon(times)
    elongation = sun_coords.separation(moon_coords)

    full_moon_idxs = argrelextrema(elongation, np.greater)
    full_moon_mjds = mjds[full_moon_idxs]
    full_moon_nights = full_moon_mjds - start_mjd
    
    return full_moon_nights

def read_downtime_lunations(file_name):
    with sqlite3.connect(file_name) as con:
        ref_downtime = pd.read_sql_query('SELECT * FROM Downtime', con)
    downtime_central_nights = ref_downtime['night'] + 0.5*ref_downtime['duration']
    lunations = np.round(downtime_central_nights.values/SYNODIC_MONTH_DAYS).astype(int)
    return lunations

def save_downtimes(downtimes, file_name):
    # Need to create table explicitly because the table
    # pandas creates automatically does not support setting
    # the primary key
    with sqlite3.connect(file_name) as con:
        cur = con.cursor()
        cur.execute("CREATE TABLE Downtime(night INTEGER PRIMARY KEY,duration INTEGER,activity TEXT);")
        
    with sqlite3.connect(file_name) as con:
        downtimes.to_sql('Downtime', con, if_exists='append', index=False)


def main():
    full_moon_nights = compute_full_moon_nights(LOCATION, START_MJD)
    downtime_lunations = read_downtime_lunations(REF_DOWNTIME_DB)
    nights = np.floor(full_moon_nights[downtime_lunations] - DOWNTIME_DURATION/2).astype(int)
    downtime = pd.DataFrame({'night': nights, 'duration': DOWNTIME_DURATION, 'activity': 'maintenance'})
    save_downtimes(downtime, 'twoweek_downtimes_a.db')
    
        
if __name__ == '__main__':
    main()
