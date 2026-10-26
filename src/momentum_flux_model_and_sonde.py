
import iris
import warnings

from iris.warnings import (IrisLoadWarning, IrisCfNonSpanningVarWarning, IrisCfMissingVarWarning)

warnings.filterwarnings("ignore", category=IrisLoadWarning)
warnings.filterwarnings("ignore", category=IrisCfNonSpanningVarWarning)
warnings.filterwarnings("ignore", category=IrisCfMissingVarWarning)

import numpy as np
import sys
import subprocess
import matplotlib.pyplot as plt
import scipy
from scipy.signal import savgol_filter
from scipy.interpolate import RegularGridInterpolator
from iris.analysis import trajectory
from iris.analysis.cartography import rotate_pole
from pathlib import Path
sys.path.append("/home/users/jmtd42/TEAMx/TEAMxRadiosondes_Public/")

import momentum_flux_funcs as tx

orog_data_path_base = "/home/users/jmtd42/TEAMx/TEAMxRadiosondes_Public/src/amrit_exploration/raso_through_model_output/"

# subprocess.run(["ls", "-l", data_path_base])

# Paths for the sonde and model data
# Pick one of these by selecting the index in the paths below. 
# Or, just provide your own path, sonde IDs, and model run 


# Output file where the orography lives
data_path_orog = "orography_u-eb392_run1_0p3km.nc"

# Original model run filepath. Used for orography inside function to calculate vertical levels
original_model_data_path = "/gws/ssde/j25b/team_x/sidmba/u-eb392_run1/0p3km_v3/w_ALPS_0p3km_ERA5_BAS_MetUM_v3_2s_20250224T1200Z.nc"

sonde_ids = []

# data_path_base = "/gws/ssde/j25b/team_x/TEAMx_FLOW/tbanyard/data/radiosondes/sterzing/202502"
data_path_base = "/home/users/jmtd42/TEAMx/TEAMxRadiosondes_Public/src/amrit_exploration/momentum_fluxes/iop_17_18/iop_17_18_raw_sonde_data" 

data_dir_sonde = data_path_base
data_dir_model = data_path_base

sonde_list = sorted(Path(data_dir_sonde).glob("*.nc"))
model_list = sorted(Path(data_dir_model).glob("*.nc"))

model_and_sonde = True
save_plot = False

for i in range(0, len(sonde_list)):
    try:
        # Paths to output the flux arrays
        flux_output_path = f"{sonde_list[i].stem}_{i}_Momentum_Flux_model"

        print(40*"###")
        print(f"PROCESSING FILE {i}: {sonde_list[i].name}")
        
        # Import in the cubes
        plot_title_affix = f"{i}"
        cubelist_s_path = f"{sonde_list[i]}"
        cubelist_m_path = f"{model_list[i]}"
        cubelist_mo_path = f"{original_model_data_path}"
        orog_cube_path = f"{orog_data_path_base}/{data_path_orog}" 

        Flux_u_s = tx.flux_from_sonde(cubelist_s_path, "zonal")
        Flux_v_s = tx.flux_from_sonde(cubelist_s_path, "meridional")   
        Flux_u_m = tx.flux_from_sonde(cubelist_m_path, "zonal")       
        Flux_v_m = tx.flux_from_sonde(cubelist_m_path, "meridional")
        alts, alts_mask = tx.nonans_cubedata(cubelist_s_path, "altitude")

        # altitude in km
        alt_s_km = alts * 1e-3

        print(f"Flux calculated, now saving array")
        # Save the flux arrays as numpy arrays
        print(flux_output_path)
        if model_and_sonde == True:
            print(f"Saving both model and sonde fluxes at {flux_output_path}")
            np.savez(f"{flux_output_path}", 
                    Flux_u_s = Flux_u_s.filled(np.nan),
                    Flux_v_s = Flux_v_s.filled(np.nan),
                    Flux_u_m = Flux_u_m.filled(np.nan),
                    Flux_v_m = Flux_v_m.filled(np.nan),
                    alt_s_km = alt_s_km.filled(np.nan))

        elif model_and_sonde == False:
            print(f"Saving Sonde data at {flux_output_path}")
            np.savez(f"{flux_output_path}", 
                Flux_u_s = Flux_u_s.filled(np.nan),
                Flux_v_s = Flux_v_s.filled(np.nan),
                alt_s_km = alt_s_km.filled(np.nan))

        else:
            print(f"Saving both model and sonde")
            np.savez(f"{flux_output_path}", 
                    Flux_u_s = Flux_u_s.filled(np.nan),
                    Flux_v_s = Flux_v_s.filled(np.nan),
                    Flux_u_m = Flux_u_m.filled(np.nan),
                    Flux_v_m = Flux_v_m.filled(np.nan),
                    alt_s_km = alt_s_km.filled(np.nan))           
        print(f"Array Saved, now moving on to next file")
        print(40*"##")

    except ValueError as no_data_found:
        if "array of sample points is empty" in str(no_data_found):
            print(f"Data missing in file. Skipping: ") 
            continue
        raise
