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

def calc_vframe(gbt, *args):
    '''
    Calculate the VFRAME for a given integration

    Parameters
    ----------
        gbt : astropy.coordinates.EarthLocation
            An astropy object representing the GBT on Earth
        t_mjd : float
            the time(s) for which to calculate VFRAME
        RA  : float
            the RA (deg J2000) of the telescope pointing 
        Dec : float
            the Dec (deg J2000) of the telescope pointing 
        veldef : str
            the reference frame for which to calculate VFRAME. Must be *-BAR, *-HEL, or *-LSR

    Returns
    ----------
        frame_vels : float or list 
            the frame velocities (m/s) for each t_mjd
    '''
    # Sanitize the arguments
    t_mjd, RA, Dec, veldef = sanitize_inputs(*args, sanitype='series')

    t = Time(t_mjd, format='mjd', scale='utc')                                                                      # [3.1]
    source_sph = SkyCoord(RA, Dec, frame='icrs', unit='deg')                                                        # [3.2]
    source_cart = source_sph.icrs.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation)   # [3.3]
    source_cart = np.array([-sc.xyz for sc in source_cart]).T

    indx_bar = veldef.str.endswith('BAR')
    indx_hel = veldef.str.endswith('HEL')
    indx_lsr = veldef.str.endswith('LSR')
    indx_top = veldef.str.endswith('TOP')
    frame_vels = pd.Series(np.ones(len(t_mjd)) * 1E50)

    # [3.4] BARYCENTRIC
    if sum(indx_bar) > 0:
        epos, evel = get_body_barycentric_posvel('earth', t[indx_bar])           # [3.4.1]
        gbt_pos, gbt_vel = gbt.get_gcrs_posvel(t[indx_bar])                      # [3.4.2]
        bvel_bg = (evel + gbt_vel) * 1000 / (24*60*60)  # [3.4.3]
        bvel_bg = np.array([(bi * u.d * u.m / u.s / u.km).xyz for bi in bvel_bg])
        frame_vels[indx_bar] = np.matmul(bvel_bg, source_cart)                         # [3.4.4]

    # [3.5] HELIOCENTRIC
    if sum(indx_hel) > 0:
        frame_vels[indx_hel] = -source_sph.radial_velocity_correction('heliocentric', obstime=t[indx_hel], location=gbt) * 1000 / (24*60*60) * u.d * u.m / u.s / u.km # [3.5.4]

    # [3.6] LSRK
    if sum(indx_lsr) > 0:
        epos, evel = get_body_barycentric_posvel('earth', t[indx_lsr])                                                                  # [3.6.1]
        gbt_pos, gbt_vel = gbt.get_gcrs_posvel(t[indx_lsr])                                                                             # [3.6.1]
        lsrk_pos = SkyCoord(ra="18:00:00", dec="+30:00:00", unit=(u.hourangle, u.deg), frame='icrs', radial_velocity=20*u.km/u.s, equinox="J1900")       # [3.6.2]
        bvel_bg = np.add(np.add(evel, gbt_vel), lsrk_pos.velocity) *1000 / (24*60*60)                                  # [3.6.3]
        bvel_bg = np.array([(bi * u.d * u.m / u.s / u.km).xyz for bi in bvel_bg])
        frame_vels[indx_lsr] = np.matmul(bvel_bg, source_cart)                                                                                # [3.6.4]

    # TOPOCENTRIC (to verify that these don't get values changed if passed)
    if sum(indx_top) > 0:
        frame_vels[indx_top] = np.zeros(len(indx_top))
    return frame_vels

def calc_lo1freq(*args):
    '''
    Calculate the LO1FREQ

    Parameters
    ----------
        f0 : float 
            the rest frequency of the source (Hz)
        rvsys : float
            the RVSYS (m/s)
        lomult : int
            the LOMULT
        looffset : float
            the LOOFFSET
        iffreq : float
            the IFFREQ (Hz)
        sideband : str
            the sideband ("UPPER" or "LOWER")
    
    Returns
    ----------
        lo1freq : float 
            the LO1FREQ
    '''
    # Sanitize the arguments
    f0, rvsys, lomult, looffset, iffreq, sideband = sanitize_inputs(*args)
    # Get the sky frequency corresponding to RVSYS
    sky_freq = v2f(f0, rvsys)
    # Turn the sky frequency into an LO1FREQ
    lo1freq = sky2lo(sky_freq, lomult, looffset, iffreq, sideband)
    return lo1freq

def calc_f_offset(*args):
    '''
    Compute the final Doppler shift of the data

    Parameters
    ----------
        lo1_freq_og : float
            the original LO1FREQ
        lo1_freq_new : float
            the new LO1FREQ
        lo_mult : int
            the LOMULT
        lo_offset : int
            the LOOFFSET
        iffreq : float
            the IFFREQ
        sideband : str
            the sideband ("UPPER" or "LOWER")
    
    Returns
    ----------
        f_sky_new : float 
            the Doppler-shifted emission frequency (using the selected velocity definition) in the selected rest frame (Hz)
        d_f_sky : float
            the difference between the new and original sky frequencies
    '''
    # Sanitize the arguments
    lo1_freq_og, lo1_freq_new, lo_mult, lo_offset, iffreq, sideband = sanitize_inputs(*args)

    f_sky_og = lo2sky(lo1_freq_og, lo_mult, lo_offset, iffreq, sideband)
    f_sky_new = lo2sky(lo1_freq_new, lo_mult, lo_offset, iffreq, sideband)
    d_f_sky = np.subtract(f_sky_new, f_sky_og)
    return f_sky_new, d_f_sky

def calc_channel_offset(*args):
    '''
    Compute the final Doppler shift of the data

    Parameters
    ----------
        f_offset : float 
            the Doppler shifted emission frequency (using the selected velocity definition) in the selected rest frame (Hz)
        channel_bw  : float
            the channel bandwidth (in Hz)
        sideband : str
            the sideband ("UPPER" or "LOWER")
    
    Returns
    ----------
        channel_offset : int 
            the number of VEGAS channels to shift by, rounded to the nearest integer
    '''
    # Sanitize the arguments
    f_offset, channel_bw, sideband = sanitize_inputs(*args)

    indx_lsb = sideband == "LOWER"
    indx_usb = sideband == "UPPER"
    channel_offset = pd.Series(np.ones(len(sideband)) * 1E50)
    if sum(indx_lsb) > 0: 
        channel_offset[indx_lsb] = np.rint(np.divide(np.array(-f_offset[indx_lsb]), np.array(channel_bw[indx_lsb])))
    if sum(indx_usb) > 0:
        channel_offset[indx_usb] = np.rint(np.divide(np.array(f_offset[indx_usb]), np.array(channel_bw[indx_usb])))
    return channel_offset

def sky2lo(*args):
    '''
    Convert a sky frequency to an LO1FREQ

    Parameters
    ----------
        f_sky : float 
            the sky frequency
        lomult : float
            the LOMULT
        looffset : float
            the LOOFFSET
        iffreq : float
            the IFFREQ
        sideband : str
            the sideband ("lower" or "upper")
    
    Returns
    ----------
        lo1freq : float 
            the LO1 frequency corresponding to the sky frequency
    '''
    # Sanitize the arguments
    f_sky, lomult, looffset, iffreq, sideband = sanitize_inputs(*args)

    indx_lsb = sideband == "LOWER"
    indx_usb = sideband == "UPPER"
    lo1freq = pd.Series(np.ones(len(sideband)) * 1E50)
    if sum(indx_lsb) > 0:
        lo1freq[indx_lsb] = np.subtract(np.divide(np.add(iffreq[indx_lsb], f_sky[indx_lsb]), lomult[indx_lsb]), looffset[indx_lsb])
    if sum(indx_usb) > 0:
        lo1freq[indx_usb] = np.subtract(np.divide(np.subtract(f_sky[indx_usb], iffreq[indx_usb]), lomult[indx_usb]), looffset[indx_usb])
    return lo1freq

def lo2sky(*args):
    '''
    Convert an LO1FREQ to a sky frequency

    Parameters
    ----------
        lo1freq : float 
            the LO1FREQ
        lomult : float
            the LOMULT
        looffset : float
            the LOOFFSET
        iffreq : float
            the IFFREQ
        sideband : str
            the sideband ("UPPER" or "LOWER")
    
    Returns
    ----------
        f_sky : float 
            the sky frequency
    '''
    # Sanitize the arguments
    lo1freq, lomult, looffset, iffreq, sideband = sanitize_inputs(*args)

    indx_lsb = sideband == "LOWER"
    indx_usb = sideband == "UPPER"
    f_sky = pd.Series(np.ones(len(sideband)) * 1E50)
    f_sky[indx_lsb] = np.subtract(np.multiply(np.add(lo1freq[indx_lsb], looffset[indx_lsb]), lomult[indx_lsb]), iffreq[indx_lsb])
    f_sky[indx_usb] = np.add(np.multiply(np.add(lo1freq[indx_usb], looffset[indx_usb]), lomult[indx_usb]), iffreq[indx_usb])
    return f_sky

def calc_rvsys(*args):
    '''
    Calculate the relativistic velocity of the source (RVSYS)
    See page 6: https://www.gb.nrao.edu/GBT/MC/doc/dataproc/gbtLOFits/gbtLOFits.pdf

    Parameters
    ----------
        vframe_old : float 
            The originally-recorded radial velocity of the reference frame (m/s)
        vframe_new : float 
            The corrected radial velocity of the reference frame (m/s)
        rvsys_og : float
            the original RVSYS value (m/s)
    
    Returns
    ----------
        rvsys : float 
            The RVSYS value (m/s)
    '''
    # Sanitize the arguments
    vframe_old, vframe_new, rvsys_og = sanitize_inputs(*args)
    # Get the original source velocity
    v_source = add_relativistic_velocities(rvsys_og, -vframe_old)
    # Add the velocities together
    rvsys = add_relativistic_velocities(v_source, vframe_new)
    return rvsys

def add_relativistic_velocities(*args):
    '''
    Adds two relativistic velocities together using the velocity-addition formula
    .. math:: \frac{v_{1} + v_{2}}{1 + \frac{v_{1}v_{2}}{c^2}}

    Parameters
    ----------
        v1 : float 
            The first velocity (m/s)
        v2 : float 
            The second velocity (m/s)
    
    Returns
    ----------
        v_total : float 
            The combined velocity (m/s)
    '''
    # Sanitize the arguments
    v1, v2 = sanitize_inputs(*args)
    c = const.c.value
    num = np.add(v1, v2)
    denom = np.add(1, np.multiply(v1, v2)/(c**2))
    v_total = np.divide(num, denom)
    return v_total

def v2f(*args):
    '''
    Calculate the frequency corresponding to a relativistic Doppler shift

    Parameters
    ----------
        f0 : float 
            The rest frequency of the source (Hz)
        v : float 
            The velocity of the source (m/s)
    
    Returns
    ----------
        f : float 
            The Doppler-shifted value of f0 (Hz)
    '''
    # Sanitize the arguments
    f0, v = sanitize_inputs(*args)

    f_num = np.subtract(1, np.divide(v,const.c.value))
    f_denom = np.add(1, np.divide(v,const.c.value))
    f = np.multiply(np.sqrt(np.divide(f_num, f_denom)), pd.Series(f0))
    return f

def sanitize_inputs(*args, sanitype='values'):
    '''
    Sanitizes inputs to make sure they're the right ty

    Parameters
    ----------
        *args : tuple 
            The original arguments
        sanitype : str 
            The output type ("values" or "series")
    
    Returns
    ----------
        sanitized_args : tuple
            The sanitized arguments
    '''
    arg_dict = {}
    for ai, a in enumerate(args):
        if type(a) == pd.DataFrame:
            for aj in a.columns:
                arg_dict[aj] = a[aj]
        else:
            arg_dict[f'arg_{ai}'] = a
    arg_df = pd.DataFrame(arg_dict)
    if sanitype == 'series':
        sanitized_args = tuple([arg_df[ac] for ac in arg_df.columns])
    elif sanitype == 'values':
        sanitized_args = tuple([arg_df[ac].values for ac in arg_df.columns])
    return sanitized_args
