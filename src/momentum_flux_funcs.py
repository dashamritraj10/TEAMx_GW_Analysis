# iris imports and warning suppression
import iris
from iris.analysis import trajectory
from iris.analysis.cartography import rotate_pole

import warnings
from iris.warnings import (IrisLoadWarning, IrisCfNonSpanningVarWarning, IrisCfMissingVarWarning)

warnings.filterwarnings("ignore", category=IrisLoadWarning)
warnings.filterwarnings("ignore", category=IrisCfNonSpanningVarWarning)
warnings.filterwarnings("ignore", category=IrisCfMissingVarWarning)

# system level imports
import os
import sys
import subprocess

# Matplotlib imports
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import AutoLocator, MaxNLocator, AutoMinorLocator

# Scipy and numpy imports
import scipy
from scipy.signal import savgol_filter
from scipy.interpolate import RegularGridInterpolator
import numpy as np

# Path adjustments
from pathlib import Path
sys.path.append("/home/users/jmtd42/TEAMx/TEAMxRadiosondes_Public/")


def fill_nans(array):
    """Return an array with the nans linearly interpolated, requires input at unmasked array"""
    mask = np.isnan(array) # create boolean mask of where nans exist
    arr_nonans = array.copy() # create a copy so as to not change original array

    # np.flatnonzero returns the coords of where argument condition is met
    # this is saying that 
    # at the location where the nans exist (np.flatnonzero(mask)), interpolate the values of 
    # the places where nans dont exist (np.flatnonzero(~mask)), with these valid  values given by array[~mask]
    # this will return an array of fill values for all the points array[mask]. Replace arr_nonans[mask]
    arr_nonans[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), array[~mask])

    return arr_nonans, mask


def vert_levels_for_sonde(sonde_data_path, model_sonde_data_path, model_data_path, model_orog_path):
    """
    Determine the model vertical-level heights crossed by a radiosonde along its trajectory.

    The radiosonde latitude and longitude coordinates are rotated into the model grid
    coordinate system and used to interpolate the model orography along the sonde path.
    The geometric height of each model vertical level is then calculated from the
    reference ellipsoid height and fraction of orographic height. These level heights
    are compared with the radiosonde altitude profile to identify where the sonde
    crosses each model level.

    Parameters
    ----------
    sonde_data_path : str
        Path to the radiosonde data file.
    model_sonde_data_path : str
        Path to the model-sonde data file.
    model_data_path : str
        Path to the model data file containing the vertical-coordinate information.
    model_orog_path : str
        Path to the model orography data file.

    Returns
    -------
    cross_heights : list
        List containing the height at which the radiosonde crosses each model vertical
        level. If a model level is not crossed by the radiosonde, the corresponding
        element is ``np.nan``.
    """
    # Import in some data to see what it looks like
    cubelist_s = iris.load(f"{sonde_data_path}")
    cubelist_m = iris.load(f"{model_sonde_data_path}")
    cubelist_mo = iris.load(f"{model_data_path}")
    orog_cube = iris.load(f"{model_orog_path}")[0]

    # Get sonde altitudes
    alt_s_cube = cubelist_s.extract_cube("altitude")
    alt_s_raw = alt_s_cube.data
    alt_s_raw.mask = False
    alt_s_nonans, alt_s_mask = fill_nans(alt_s_raw)

    # Sample points to slice the lat at
    lat_cube = cubelist_s.extract_cube("latitude")
    lon_cube = cubelist_s.extract_cube("longitude")

    lats = lat_cube.data
    lons = lon_cube.data

    lats.mask = False
    lons.mask = False

    # Get grid_lat lons
    grid_lats = orog_cube.coord("grid_latitude").points
    grid_lons = orog_cube.coord("grid_longitude").points

    # Convert the points to cubes coord system
    coord_system = orog_cube.coord("grid_latitude").coord_system

    # Find the pole
    pol_lat = coord_system.grid_north_pole_latitude
    pol_lon = coord_system.grid_north_pole_longitude

    # Rotate the coords, and remove nans
    rot_lats_raw, rot_lons_raw = rotate_pole(lats=lats, lons=lons, pole_lat=pol_lat, pole_lon=pol_lon)
    rot_lons, lon_mask = fill_nans(rot_lons_raw[::] + 360)
    rot_lats, lat_mask = fill_nans(rot_lats_raw[::])

    # Create the sample points to slice cube at
    grid_sample_points = (grid_lats, grid_lons)
    sonde_sample_points = np.column_stack((rot_lats, rot_lons))

    orog_data_raw = orog_cube.data
    orog_data = orog_data_raw[0,0,::,::]
    orog_data.mask=False

    # Create interpolator and use it 
    orog_interpolator = RegularGridInterpolator(points=grid_sample_points, values=orog_data, method="linear")
    orog = orog_interpolator(sonde_sample_points)

    # Now to figure out how to get information about the vertical levels out
    mlevels_cube = cubelist_mo.extract_cube("model_level_number")
    surf_alt_cube = cubelist_mo.extract_cube("surface_altitude")
    w_cube = cubelist_mo.extract_cube("upward_air_velocity")

    # Extract the auxillary coord data
    orog_frac = w_cube.aux_coords[0].points
    ref_ellip_height = w_cube.aux_coords[1].points

    # Level heights
    lev_index = mlevels_cube.data

    levs_list = []
    for i,j in enumerate(lev_index):
        lev_height = ref_ellip_height[i] + orog_frac[i]*orog
        levs_list.append(lev_height)
        
    cross_heights = []

    for vert_level in levs_list:
    
        # for each time alt gives the vertical location of the sonde
        # for each time, vert_level gives the height of this vertical level at 
        # the horizontal position the sonde is at
        # the difference quantifies if the sonde is below, at, or above this height  
        diff = alt_s_nonans - vert_level

        # The sign of the difference tells us if the sonde is below, at, or above
        # the height of this vertical level at this time.
        # 
        # Then, investigate where these signs change using np.diff 
        # Where is this sign change != 0?, i.e. where does sonde cross
        # this vertical level?
        # 
        # The output is a mask, telling us the location of each time stamp 
        # that the sonde crosses a vertical level 
        # 
        
        crossing_mask = np.where(np.diff(np.sign(diff))!=0)

        # Now, using this mask pick out the height in vert_level 
        # where the crossing has occured

        crossing_point = vert_level[crossing_mask]

        if np.size(crossing_point) == 0:
            crossing_point = np.nan
        
        cross_heights.append(crossing_point)

    return cross_heights


def nonans_cubedata(cubepath, var_name):

    # Get the list of cubes from the path
    cubelist = iris.load(cubepath)

    # Extract the cube
    var_cube = cubelist.extract_cube(var_name)

    # Remove nans and masks from the data
    var_cube_raw = var_cube.data
    var_cube_raw.mask = False
    var_cube_nonans, var_cube_mask = fill_nans(var_cube_raw)

    return var_cube_nonans, var_cube_mask


def coriolis_param_for_sonde(cubepath):
    lat_nonans, lat_mask = nonans_cubedata(cubepath, "latitude")
    lats_median = np.median(lat_nonans)
    coriolis_param = coriolis_param_at_lat(lats_median)

    return coriolis_param

def coriolis_param_at_lat(lat, lat_units="degrees", g=9.81, Omega=7.2921 * 1e-5):
    # Use default params to find inertial frequency
    if lat_units == "degrees":
        lat_rads = np.pi/180 * lat
    else: lat_rads = lat

    coriolis_param = 2 * Omega * np.sin(lat_rads)

    return coriolis_param


def get_means_primes_pfit(var_array, alt_array, p_deg=4):
    var_coeffs = np.polynomial.Polynomial.fit(alt_array, var_array, p_deg)
    var_mean = var_coeffs(alt_array)

    var_prime = var_array - var_mean

    return var_prime, var_mean, var_coeffs


def plot_flux(plot_title_affix, cubelist_s_path, cubelist_m_path, 
              cubelist_mo_path, orog_cube_path, save_plot = True):

    
    cross_heights = vert_levels_for_sonde(sonde_data_path=cubelist_s_path, 
                                        model_sonde_data_path=cubelist_m_path, 
                                        model_data_path=cubelist_mo_path, 
                                        model_orog_path=orog_cube_path)
        

    Flux_u_s = flux_from_sonde(cubelist_s_path, "zonal")
    Flux_v_s = flux_from_sonde(cubelist_s_path, "meridional")   
    Flux_u_m = flux_from_sonde(cubelist_m_path, "zonal")       
    Flux_v_m = flux_from_sonde(cubelist_m_path, "meridional")
    alts_km = nonans_cubedata(cubelist_s_path, "altitude")


# Plot the Fluxes
    fig = plt.figure(figsize=(10,8), layout="constrained") 
    ax = fig.subplots(1,2).flatten()

    ylim = (0,None)

    max_flux = np.max(np.array([np.max(Flux_u_s), np.max(Flux_u_m), np.max(Flux_v_s), np.max(Flux_v_m)]))
    min_flux = np.min(np.array([np.min(Flux_u_s), np.min(Flux_u_m), np.min(Flux_v_s), np.min(Flux_v_m)]))

    xlim = (1.25*min_flux, 1.25*max_flux)

    # Plot zonal flux
    ax[0].plot(Flux_u_s, alts_km, color="C1")
    ax[0].plot(Flux_u_m, alts_km, color="k")
    ax[0].axvline(0, linestyle="--", linewidth=1, color="k", alpha=0.5)

    ax[0].set_xlim(xlim)
    ax[0].locator_params(axis="x", nbins=15)
    ax[0].set_xlabel(r"$\rho \overline{u^\prime w^\prime}$ ($Pa$)")

    ax[0].set_ylim(ylim)
    ax[0].locator_params(axis="y", nbins=15)
    ax[0].set_ylabel(r"Altitude ($km$)")

    ax[0].set_title(rf"Pseudomomentum Flux $F_u$ sonde: {plot_title_affix}")

    #Plot meridional flux
    ax[1].plot(Flux_v_s, alts_km, color="C1")
    ax[1].plot(Flux_v_m, alts_km, color="k")
    ax[1].axvline(0, linestyle="--", linewidth=1, color="k", alpha=0.5)

    ax[1].set_xlim(xlim)
    ax[1].locator_params(axis="x", nbins=15)
    ax[1].set_xlabel(r"$\rho \overline{v^\prime w^\prime}$ ($Pa$)")

    ax[1].set_ylim(ylim)
    ax[1].locator_params(axis="y", nbins=15)
    ax[1].set_ylabel(r"Altitude ($km$)")

    ax[1].set_title(rf"Pseudomomentum Flux $F_v$ sonde: {plot_title_affix}")

    for height in cross_heights:
        if np.isfinite(height):
            ax[0].axhline(height*1e-3, linewidth=0.3, alpha=1, zorder=1)
            ax[1].axhline(height*1e-3, linewidth=0.3, alpha=1, zorder=1)

    if save_plot == True:
        plt.savefig(f"Momentum_flux_sonde_id_{plot_title_affix}.png")
    else: 
        plt.close()

    return

def ave_spectral_frequency(N, f, p=5/3, g=9.81):
    q = p-1
    r = 2-p
    f_hat = f/N
    omega_bar = N*((q)/(r))*(f_hat**(q))*(1-f_hat**(r))

    return omega_bar

def inf_to_nan(array):
    array_inf_to_nan = array.copy()
    array_inf_to_nan[np.where(np.imag(array) !=0 )] = np.nan
    array_real = np.real(array_inf_to_nan)

    return array_real


def brunt_vaisala_from_sonde(cubelist_path, p_deg = 5, g=9.81):
    Theta_nonans, Theta_mask = nonans_cubedata(cubelist_path, "air_potential_temperature")
    alt_nonans, alt_mask = nonans_cubedata(cubelist_path, "altitude")

    Theta_coeffs = np.polynomial.Polynomial.fit(alt_nonans, Theta_nonans, p_deg)

    Theta_mean = Theta_coeffs(alt_nonans)
    Theta_prime = Theta_nonans - Theta_mean

    # Now need the Brunt Vaisala Frequency 
    dthetacoeffs = Theta_coeffs.deriv()
    dtheta_dz_s = dthetacoeffs(alt_nonans)
    N_sqrd = g/(Theta_mean) * dtheta_dz_s

    print(f"Brunt Vaisala Frequency Calculated")

    # To find N_s, ignore points where N_sqrd is negative
    N = np.full_like(N_sqrd, np.nan, dtype=float)
    valid = N_sqrd > 0 
    N[valid] = np.sqrt(N_sqrd[valid])

    return N


def rolling_mean(array, window_size):
    # Take a rolling mean to find average
    window = min(window_size, len(array))
    kernel = np.ones(window) / window
    array_mean = np.convolve(array, kernel, mode="same")

    return array_mean

def flux_from_sonde(cubelist_path, component="zonal"):
    if component == "zonal": var_name = "eastward_wind"
    elif component == "meridional": var_name = "northward_wind"
    else:
        print("Wind component not specified correctly")
        raise 

    #print(f"Loading {component} wind data and Temperature")
    wind_nonans, u_wind_mask = nonans_cubedata(cubelist_path, var_name)
    T_nonans, T_mask = nonans_cubedata(cubelist_path, "air_temperature")

    #print(f"Loading in Pressure data and converting into Pascals")
    pressure_nonans_hpa, pressure_mask = nonans_cubedata(cubelist_path, "air_pressure")
    pressure_nonans = pressure_nonans_hpa * 100

    #print(f"Loading in altitude and time data")
    alt_nonans, alt_mask = nonans_cubedata(cubelist_path, "altitude")
    time_s, time_mask = nonans_cubedata(cubelist_path, "Elapsed Time")

    #print(f"Data loaded, now fitting polynomials to find means and primes")
    # Now, fit a polynomial to the u and T data and find coeffs
    p_deg = 4
    wind_coeffs = np.polynomial.Polynomial.fit(alt_nonans, wind_nonans, p_deg)
    T_coeffs = np.polynomial.Polynomial.fit(alt_nonans, T_nonans, p_deg)
    pressure_coeffs = np.polynomial.Polynomial.fit(alt_nonans, pressure_nonans, p_deg)

    # Find the mean profiles using the calculated coefficients
    wind_mean = wind_coeffs(alt_nonans)
    T_mean = T_coeffs(alt_nonans)
    pressure_mean = pressure_coeffs(alt_nonans)
    
    # Find the perturbations using raw - mean
    wind_prime = wind_nonans - wind_mean
    T_prime = T_nonans - T_mean 

    print(f"Means and Perturbations for Wind, Temperature, Pressure found.")
    
    # Normalise the Temperatures to find T prime hat
    T_hat_prime = T_prime / T_mean

    #print(f"Found normalised Temperature")

    # Use an ideal EOS to find the densities
    rho = (pressure_mean) / (287.05 * T_mean)

    # Now find T prime + 90 by using hilbert transform
    T_hat_prime_90 = np.imag(scipy.signal.hilbert(T_hat_prime))

    # Before constructing u'T^hat_+90', smooth the Temperature to remove 
    # small scale stuff
    T_hat_prime_90_smth = rolling_mean(T_hat_prime_90, window_size=50)

    print(f"Hilbert Transformed Normalised Temperature")
    print(f"Shape of transformed temperature: {np.shape(T_hat_prime_90_smth)}")

    # Now construct the u, T, rho term for sonde 
    u_prime_T_hat_prime_rho = T_hat_prime_90_smth * wind_prime * rho

    # Take a rolling mean to find average
    uT_hat_mean = rolling_mean(u_prime_T_hat_prime_rho, window_size=50)

    #print(f"Flux term calculated, now calculating prefactor")

    # Find the inertial frequency f. For this need the latitude and Omega
    f_inertial = coriolis_param_for_sonde(cubelist_path)

    # Flux 
    N = brunt_vaisala_from_sonde(cubelist_path)
    N_sqrd = N**2
    omega_bar = ave_spectral_frequency(N, f_inertial)
    g = 9.81
    # Flux
    F_u = omega_bar/N_sqrd * g * uT_hat_mean

    # Take the infinities from the flux and make them into nans
    F_u_inf_to_nan = F_u.copy()
    F_u_inf_to_nan[np.where(np.imag(F_u) !=0 )] = np.nan
    F_u_real = np.real(F_u_inf_to_nan)

    test_factor = 1
    # Calculating flux in Pa
    Flux_u = F_u_real * test_factor

    return Flux_u


def add_profile_to_ax(plot_axis, F_mag_s, F_mag_m, f_u_dot_m, f_u_dot_s, alt_km, title=None):

    ax = plot_axis
    # Set up a norm for coloured profile sonde
    points_s = np.array([F_mag_s, alt_km]).T.reshape(-1,1,2)
    norm_s = TwoSlopeNorm(vmin=-1, vmax=1, vcenter=0)
    segments_s = np.concatenate([points_s[:-1], points_s[1:]], axis=1)

    # Create linecollection for sonde and add to plot
    lc_s = LineCollection(segments_s, cmap="coolwarm", norm=norm_s, linewidth=2)
    lc_s.set_array(f_u_dot_s)
    ax.add_collection(lc_s)

    # Set up a norm for coloured profile model
    points_m = np.array([-F_mag_m, alt_km]).T.reshape(-1,1,2)
    norm_m = TwoSlopeNorm(vmin=-1, vmax=1, vcenter=0)
    segments_m = np.concatenate([points_m[:-1], points_m[1:]], axis=1)

    # Create linecollection for model and add to plot
    lc_m = LineCollection(segments_m, cmap="coolwarm", norm=norm_m, linewidth=2)
    lc_m.set_array(f_u_dot_m)
    ax.add_collection(lc_m)

    # Set title and other ax params
    ax.set_title(title)
    x_max = np.max(np.array([np.max(F_mag_m), np.max(F_mag_s)]))
    ax.set_xlim(-1.25*x_max, 1.25*x_max)
    ax.set_ylim(0, None)

    ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # Add a dashed line for 0
    ax.axvline(0, color="k", linewidth=1, linestyle="--", alpha=0.2)

    return ax, lc_s


def add_profile_to_ax(plot_axis, F_mag_s, F_mag_m, f_u_dot_m, f_u_dot_s, alt_km, vmin=-1, vmax=1, vcenter=0 ,cmap="coolwarm" ,title=None):

    ax = plot_axis
    # Set up a norm for coloured profile sonde
    points_s = np.array([F_mag_s, alt_km]).T.reshape(-1,1,2)
    norm_s = TwoSlopeNorm(vmin=vmin, vmax=vmax, vcenter=vcenter)
    segments_s = np.concatenate([points_s[:-1], points_s[1:]], axis=1)

    # Create linecollection for sonde and add to plot
    lc_s = LineCollection(segments_s, cmap=cmap, norm=norm_s, linewidth=2)
    lc_s.set_array(f_u_dot_s)
    ax.add_collection(lc_s)

    # Set up a norm for coloured profile model
    points_m = np.array([-F_mag_m, alt_km]).T.reshape(-1,1,2)
    norm_m = TwoSlopeNorm(vmin=vmin, vmax=vmax, vcenter=vcenter)
    segments_m = np.concatenate([points_m[:-1], points_m[1:]], axis=1)

    # Create linecollection for model and add to plot
    lc_m = LineCollection(segments_m, cmap=cmap, norm=norm_m, linewidth=2)
    lc_m.set_array(f_u_dot_m)
    ax.add_collection(lc_m)

    # Set title and other ax params
    ax.set_title(title)
    x_max = np.max(np.array([np.max(F_mag_m), np.max(F_mag_s)]))
    ax.set_xlim(-1.25*x_max, 1.25*x_max)
    ax.set_ylim(0, None)

    ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # Add a dashed line for 0
    ax.axvline(0, color="k", linewidth=1, linestyle="--", alpha=0.2)

    return ax, lc_s

def return_barbs(plot_axis, u_array, v_array, alt_array, barb_count=25, barb_length=5, 
                 color="black", position=0.025, norm_mag=10):
    ax = plot_axis
    barb_locs = np.linspace(0, np.size(alt_array)-1, barb_count, dtype=int)
    s = barb_locs
    speed = np.hypot(u_array, v_array)
    u_array_hat, v_array_hat = norm_mag*u_array/speed, norm_mag*v_array/speed
    barb_u = u_array_hat[s]
    barb_v = v_array_hat[s]
    barb_alts = alt_array[s]

    ax.barbs(np.full_like(barb_alts, position), barb_alts, barb_u, barb_v, transform=ax.get_yaxis_transform(),
             pivot="tip", color=color, length=barb_length)

    return ax


def flux_wind_dot(f_u, f_v, u, v):
    F = np.column_stack((f_u, f_v))
    U = np.column_stack((u, v))

    F_mag = np.linalg.norm(F, axis=1)
    U_mag = np.linalg.norm(U, axis=1)

    f_u_dot = np.sum(F * U, axis=1) / (F_mag * U_mag)

    return f_u_dot