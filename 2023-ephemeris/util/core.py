import numpy as np
import pandas as pd
import astropy.constants as const
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, solar_system_ephemeris, get_body_barycentric_posvel, CartesianRepresentation, UnitSphericalRepresentation

# Use the JPL DE430 ephemeris 
solar_system_ephemeris.set('jpl') 

def create_obj_gbt():
    '''
    Create an astropy EarthLocation for the GBT using the same values the regular system does
    Using values from page 3: https://www.gb.nrao.edu/GBT/MC/doc/dataproc/gbtLOFits/gbtLOFits.pdf
    But could also see line 69 of SolSys/astrtime.cc for the "true" values...
        latitude    = 38d 25m 59.265s N
        longitude   = 79d 50m 23.419s W
        height      = 854.83 m 

    Returns
    ----------
        gbt : EarthLocation 
            astropy EarthLocation for the GBT
    '''
    gbt_lat     =  38.433129 * u.deg    # [1.1]
    gbt_lon     = -79.839839 * u.deg    # [1.2]
    gbt_height  =  854.83   * u.m       # [1.3]
    gbt         = EarthLocation.from_geodetic(lon=gbt_lon, lat=gbt_lat, height=gbt_height)
    return gbt

def calc_vframe(gbt, t_mjd, RA, Dec, ref_frame):
    '''
    Calculate the VFRAME for a given integration

    Parameters
    ----------
        gbt : astropy.coordinates.EarthLocation
            An astropy object representing the GBT on Earth
        t_mjd : float or list
            the time(s) for which to calculate VFRAME
        RA  : float or list 
            the RA (deg J2000) of the telescope pointing 
        Dec : float or list 
            the Dec (deg J2000) of the telescope pointing 
        ref_frame : str
            the reference frame for which to calculate VFRAME. Must be *-BAR, *-HEL, or *-LSR

    Returns
    ----------
        frame_vels : float or list 
            the frame velocities (m/s) for each t_mjd
    '''
    t = Time(t_mjd, format='mjd', scale='utc')                                                                      # [3.1]
    source_sph = SkyCoord(RA, Dec, frame='icrs', unit='deg')                                                        # [3.2]
    source_cart = source_sph.icrs.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation)   # [3.3]
    
    # [3.4] 
    if ref_frame[-3:] == 'BAR':
        epos, evel = get_body_barycentric_posvel('earth', t)                    # [3.4.1]
        gbt_pos, gbt_vel = gbt.get_gcrs_posvel(t)                               # [3.4.2]
        bvel_bg = (evel + gbt_vel) *1000 / (24*60*60) * u.d * u.m / u.s / u.km  # [3.4.3]
        frame_vels = bvel_bg.dot(-source_cart)                                  # [3.4.4]
    # [3.5]
    elif ref_frame[-3:] == 'HEL':
        frame_vels = -source_sph.radial_velocity_correction('heliocentric', obstime=t, location=gbt) * 1000 / (24*60*60) * u.d * u.m / u.s / u.km # [3.5.4]
    # [3.6]
    elif ref_frame[-3:] == 'LSR':
        epos, evel = get_body_barycentric_posvel('earth', t)                                                                            # [3.6.1]
        gbt_pos, gbt_vel = gbt.get_gcrs_posvel(t)                                                                                       # [3.6.1]
        lsrk_pos = SkyCoord(ra="18:03:50.24", dec="+30:00:16.8", unit=(u.hourangle, u.deg), frame='icrs', radial_velocity=20*u.km/u.s)  # [3.6.2]
        bvel_bg = (evel + gbt_vel + lsrk_pos.velocity) *1000 / (24*60*60) * u.d * u.m / u.s / u.km                                      # [3.6.3]
        frame_vels = bvel_bg.dot(-source_cart)                                                                                          # [3.6.4]
    else:
        print("Imvalid frame. Must be *-BAR, *-HEL, or *-LSR")
    return frame_vels.value

def calc_vframe_offset(vframe_og, vframe_new):
    '''
    Calculate how much VFRAME needs to shift based on true value and what was applied

    Parameters
    ----------
        vframe_og : float or list 
            the VFRAME from the original observation 
        vframe_new : float or list
            the calculated true VFRAME

    Returns
    ----------
        vframe_offset : numpy array
            the VFRAME offsets to use to shift VEGAS data
    '''
    vframe_offset = np.array(vframe_new) - vframe_og    # [4.1]
    return vframe_offset

def calc_f_offset(f_0, vframe, lo1_freq_og, lo_mult, iffreq, sideband, formula='rel'):
    '''
    Compute the final Doppler shift of the data

    Parameters
    ----------
        f_0 : float 
            the rest frequency for which to shift
        vframe  : float
            the new velocity of the reference frame
    
    Returns
    ----------
        f_offset : float 
            the Doppler-shifted emission frequency (using the selected velocity definition) in the selected rest frame (Hz)
    '''
    c = const.c.value
    f_offset_og = lo2sky(lo1_freq_og, lo_mult, iffreq, sideband)
    v_over_c = np.divide(vframe, c)
    # [5.2.1]
    if formula == 'rel':
        f_new = f_0 * np.divide(np.sqrt(1 - v_over_c**2), 1 + v_over_c)
    # [5.2.2]
    elif formula == 'opt':
        f_new = f_0 / (1 + v_over_c)
    # [5.2.3]
    elif formula == 'rad':
        f_new = f_0 * (1 - v_over_c)
    f_offset = f_new - f_offset_og
    return f_new, f_offset


def calc_channel_offset(f_offset, channel_bw, sideband):
    '''
    Compute the final Doppler shift of the data

    Parameters
    ----------
        f_offset : float 
            the Doppler shifted emission frequency (using the selected velocity definition) in the selected rest frame (Hz)
        channel_bw  : float
            the channel bandwidth (in Hz)
    
    Returns
    ----------
        channel_offset : int 
            the number of VEGAS channels to shift by, rounded to the nearest integer
    '''
    if sideband == "LOWER":
        channel_offset = int(np.round(-f_offset/channel_bw))
    else:
        channel_offset = int(np.round(f_offset/channel_bw))
    return channel_offset

def sky2lo(f_sky, lomult, looffset, iffreq, sideband):
    '''
    Convert LO1FREQ to a sky frequency

    Parameters
    ----------
        f_sky : float 
            the sky frequency
        lomult : np.array (3,)
            the LO multiplier
        iffreq : float
            the IF frequency
        sideband : str
            the sideband ("lower" or "upper")
    
    Returns
    ----------
        lo1freq : float 
            the LO1 frequency corresponding to the sky frequency
    '''
    if sideband == "LOWER":
        lo1freq = np.divide(np.add(iffreq, f_sky), lomult) - looffset
    elif sideband == "UPPER":
        lo1freq = np.divide(np.subtract(f_sky, iffreq), lomult) - looffset
    else:
        print(f"sky2lo: {sideband} is not a valid sideband. Value must be \"LOWER\" or \"UPPER\"")
        lo1freq = np.nan
    return lo1freq

def lo2sky(lo1freq, lomult, iffreq, sideband):
    '''
    Convert LO1FREQ to a sky frequency

    Parameters
    ----------
        lo1freq : float 
            the LO1 frequency
        lomult : np.array (3,)
            the LO multiplier
        iffreq : float
            the IF frequency
        sideband : str
            the sideband ("lower" or "upper")
    
    Returns
    ----------
        f_sky : float 
            the sky frequency
    '''
    if sideband == "LOWER":
        f_sky = np.subtract(np.multiply(lo1freq, lomult), iffreq)
    elif sideband == "UPPER":
        f_sky = np.add(np.multiply(lo1freq, lomult), iffreq)
    else:
        print(f"sky2lo: {sideband} is not a valid sideband. Value must be \"LOWER\" or \"UPPER\"")
        f_sky = np.nan
    return f_sky

def calc_rvsys(v, vframe, veldef):
    '''
    Calculate RVSYS
    See page 6: https://www.gb.nrao.edu/GBT/MC/doc/dataproc/gbtLOFits/gbtLOFits.pdf

    Parameters
    ----------
        f_sky : float 
            the sky frequency
        lomult : np.array (3,)
            the LO multiplier
        iffreq : float
            the IF frequency
        sideband : str
            the sideband ("lower" or "upper")
    
    Returns
    ----------
        lo1freq : float 
            the LO1 frequency corresponding to the sky frequency
    '''
    veldef_series = pd.Series(veldef)
    opt_indx = veldef_series.str.startswith('VOPT')
    rad_indx = veldef_series.str.startswith('VRAD')
    rel_indx = veldef_series.str.startswith('VREL')

    v_corr = np.ones(len(v))
    v_corr[opt_indx] = np.divide(v[opt_indx], const.c.value) + 1
    v_corr[rad_indx] = 1 - np.divide(v[rad_indx], const.c.value)
    souvel_num = 1 + v_corr
    souvel_denom = 1 + np.multiply(v_corr, v_corr)
    souvel = np.multiply(v, np.divide(souvel_num, souvel_denom))

    r_num = np.add(souvel, vframe)
    r_denom = 1 + np.multiply(souvel, vframe)/(const.c.value**2)
    rvsys = np.divide(r_num, r_denom)
    return rvsys

