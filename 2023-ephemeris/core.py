import numpy as np
import astropy.constants as const
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, solar_system_ephemeris, get_body_barycentric_posvel, CartesianRepresentation, UnitSphericalRepresentation
from scipy.interpolate import interp1d
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
    gbt_lat     =  38.433129 * u.deg
    gbt_lon     = -79.839839 * u.deg
    gbt_height  =  8254.83   * u.m
    gbt         = EarthLocation.from_geodetic(lon=gbt_lon, lat=gbt_lat, height=gbt_height)
    return gbt

def calc_vframe(t_mjd, RA, Dec, ref_frame):
    '''
    Calculate the VFRAME for a given integration

    Parameters
    ----------
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
    t = Time(t_mjd, format='mjd', scale='utc')
    gbt = create_obj_gbt()
    source_sph = SkyCoord(RA, Dec, frame='icrs', unit='deg')#, radial_velocity=df_bar_b['velocity'] * u.km / u.s)
    source_cart = source_sph.icrs.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation)
    
    if ref_frame[-3:] == 'BAR':
        epos, evel = get_body_barycentric_posvel('earth', t)
        gbt_pos, gbt_vel = gbt.get_gcrs_posvel(t)
        bvel_bg = (evel + gbt_vel) *1000 / (24*60*60) * u.d * u.m / u.s / u.km
        frame_vels = bvel_bg.dot(-source_cart)
    elif ref_frame[-3:] == 'HEL':
        frame_vels = -source_sph.radial_velocity_correction('heliocentric', obstime=t, location=gbt) * 1000 / (24*60*60) * u.d * u.m / u.s / u.km
    elif ref_frame[-3:] == 'LSR':
        lsrk_pos = SkyCoord(ra="18:03:50.24", dec="+30:00:16.8", unit=(u.hourangle, u.deg), frame='icrs', radial_velocity=20*u.km/u.s)
        epos, evel = get_body_barycentric_posvel('earth', t)
        gbt_pos, gbt_vel = gbt.get_gcrs_posvel(t)
        bvel_bg = (evel + gbt_vel + lsrk_pos.velocity) *1000 / (24*60*60) * u.d * u.m / u.s / u.km
        frame_vels = bvel_bg.dot(-source_cart)
    else:
        print("Imvalid frame. Must be *-BAR, *-HEL, or *-LSR")
    return frame_vels.value

def calc_vframe_offset(DMJD_og, DMJD_new, vframe_og, vframe_new):
    '''
    Calculate how much VFRAME needs to shift based on true value and what was applied

    Parameters
    ----------
        DMJD_og : float or list
            the DMJD of the original VFRAME(s)
        DMJD_new  : float or list 
            the DMJD of each VEGAS integration 
        vframe_og : float or list 
            the VFRAME from the original observation 
        vframe_new : float or list
            the calculated true VFRAME

    Returns
    ----------
        vframe_offset : numpy array
            the VFRAME offsets to use to shift VEGAS data
    '''
    vframe_offset = np.array(vframe_new) - vframe_og
    return vframe_offset

def calc_f_offset(f_0, vframe_offset, formula='rel'):
    '''
    Compute the final Doppler shift of the data

    Parameters
    ----------
        f_0 : float 
            the rest frequency for which to shift
        vframe_offset  : float
            the leftover velocity of the reference frame that needs to be added
    
    Returns
    ----------
        f_offset : float 
            the Doppler shifted emission frequency (using the selected velocity definition) in the selected rest frame (Hz)
    '''
    #f_0 = f_0 * u.Hz
    c = const.c.value
    #vframe_offset = vframe_offset * u.m / u.s
    if formula == 'rel':
        # swapped minus and plus
        f_offset = f_0 * np.sqrt(np.divide((c-vframe_offset),(c+vframe_offset))) - f_0
    elif formula == 'opt':
        f_offset = f_0 / (1 + vframe_offset/c) - f_0
    elif formula == 'rad':
        f_offset = f_0 * (1 - vframe_offset/c) - f_0
    return -f_offset


def calc_channel_offset(f_offset, channel_bw):
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
    channel_offset = int(np.round(f_offset/channel_bw))
    return channel_offset

def sky2lo(f_sky, lomult, iffreq, vframe, sideband, formula='rel', vel=0):
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
    source_vel_f_off = calc_f_offset(f_sky, vel, formula) 
    f_sky -= source_vel_f_off
    source_vel_f_off = calc_f_offset(f_sky, vframe, 'rel') 
    f_sky -= source_vel_f_off
    if sideband == "LOWER":
        lo1freq = np.divide(np.add(iffreq, f_sky), lomult)
    elif sideband == "UPPER":
        lo1freq = np.divide(np.subtract(f_sky, iffreq), lomult)
    else:
        print(f"sky2lo: {sideband} is not a valid sideband. Value must be \"LOWER\" or \"UPPER\"")
        lo1freq = np.nan
    return lo1freq

def calc_rvsys(souvel, vframe):
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
    rvsys = np.divide(np.add(souvel, vframe), 1 + np.multiply(souvel, vframe)/(const.c.value**2))
    return rvsys

