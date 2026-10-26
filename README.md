# TEAMx Gravity Wave Analysis

Python scripts and Jupyter notebooks for processing radiosonde and numerical model data, calculating gravity-wave momentum fluxes, and comparing observational and model-derived estimates.

**Status:** Research code — documentation and organisation in progress.

## 1. Repository Structure

```text
TEAMx_GW_Analysis/
│
├── README.md
├── requirements.txt
│
├── src/
│   ├── momentum_flux.py
│   ├── momentum_flux_funcs.py
│   └── momentum_flux_model_and_sonde.py
│
├── notebooks/
│   └── *.ipynb
│
├── data/
│   ├── IOP_dates.txt
│   └── iop_17_18/
│       ├── iop_17_18_raw_sonde_data/
│       ├── iop_17_18_coarse_sonde_data/
│       ├── iop_17_18_model_sonde_data/
│       ├── iop_17_18_raw_sonde_flux_arrays/
│       ├── iop_17_18_coarse_sonde_flux_arrays/
│       ├── iop_17_18_model_sonde_flux_arrays/
│       └── orography_u-eb392_run1_0p3km.nc
│
└── results/
    └── output_plots/
```

## 2. Source Code (`src/`)

The main analysis code is contained in three Python files.

| File | Description |
|---|---|
| `momentum_flux_funcs.py` | Supporting functions for processing atmospheric profiles and calculating momentum-flux diagnostics. |
| `momentum_flux.py` | Analysis and calculation of momentum fluxes from atmospheric profiles. |
| `momentum_flux_model_and_sonde.py` | Processing and comparison of momentum-flux estimates from radiosonde and model data. |

### Code Organisation

The intended distinction between the files is:

- **`momentum_flux_funcs.py`** — Reusable functions implementing individual calculations.
- **`momentum_flux.py`** — Code performing the momentum-flux analysis.
- **`momentum_flux_model_and_sonde.py`** — Code applying the analysis to model and radiosonde datasets.

The precise import relationships, function signatures, and execution order remain to be documented.

## 3. Data Organisation (`data/`)

Data are organised by Intensive Observation Period (IOP).

Currently, the repository contains data associated with IOPs 17 and 18.

### Input Data

| Directory | Contents | Format |
|---|---|---|
| `iop_17_18_raw_sonde_data/` | Original radiosonde observations | `.nc` |
| `iop_17_18_coarse_sonde_data/` | Radiosonde data represented on a coarser vertical grid | `.nc` |
| `iop_17_18_model_sonde_data/` | Model atmospheric profiles associated with radiosonde launches | `.nc` |

Additional input files:

- `IOP_dates.txt` — List of observation periods and associated dates.
- `orography_u-eb392_run1_0p3km.nc` — Model orography data.

### Processed Data

The momentum-flux calculations generate NumPy archives containing the calculated arrays.

These are stored in:

```text
iop_17_18_raw_sonde_flux_arrays/
iop_17_18_coarse_sonde_flux_arrays/
iop_17_18_model_sonde_flux_arrays/
```

Each directory contains `.npz` files corresponding to individual radiosonde launches.

The filenames retain identifiers from the associated input datasets, allowing the processed arrays to be matched with their source profiles.

The contents and array keys of the `.npz` archives remain to be documented.

## 4. Analysis Workflow

The code broadly follows the workflow:

```text
                 INPUT DATA
                     |
          +----------+----------+
          |          |          |
       Raw sonde  Coarse sonde  Model
          |          |          |
          +----------+----------+
                     |
              Data processing
                     |
           Background profiles
                     |
           Perturbation analysis
                     |
          Momentum-flux calculation
                     |
              Saved .npz arrays
                     |
          Comparison and plotting
                     |
               Output figures
```

### Step 1 — Load Atmospheric Profiles

The analysis uses NetCDF files containing radiosonde measurements and corresponding model profiles.

The data include atmospheric quantities such as:

- Zonal and meridional wind components.
- Temperature and potential temperature.
- Pressure.
- Altitude.
- Geographical coordinates.

### Step 2 — Calculate Background Profiles

Background atmospheric profiles are estimated using polynomial fitting.

Perturbations are calculated by subtracting the fitted background from the original profiles.

### Step 3 — Calculate Momentum Fluxes

The momentum-flux calculations use atmospheric perturbations and associated diagnostics, including the Brunt–Väisälä frequency and Coriolis parameter.

The observational analysis includes calculations involving Hilbert-transformed temperature perturbations.

### Step 4 — Save Calculated Arrays

Calculated quantities are saved as `.npz` archives.

The outputs are separated according to whether they originate from raw radiosonde observations, coarse-resolution observations, or model profiles.

### Step 5 — Compare and Visualise Results

The Jupyter notebooks load and analyse atmospheric profiles and previously calculated flux arrays.

They are used to produce diagnostic plots, compare different datasets, and investigate individual aspects of the results.

**Note:** This describes the conceptual workflow. The exact execution sequence and the scripts responsible for producing each intermediate output have not yet been verified.

## 5. Jupyter Notebooks (`notebooks/`)

The notebooks are primarily used for exploratory analysis and visualisation.

| Notebook | Purpose |
|---|---|
| `coarse_sonde_vs_raw_sonde.ipynb` | Compare original and coarse-resolution radiosonde profiles. |
| `flux_s_m_comparison.ipynb` | Compare model and radiosonde momentum-flux estimates. |
| `fluxes_entire_eop.ipynb` | Analyse momentum fluxes across an extended observation period. |
| `pseudo_momentum_flux.ipynb` | Investigate pseudo-momentum-flux calculations. |
| `scorer_parameter.ipynb` | Calculate or investigate the Scorer parameter. |
| `sonde_3D_plots.ipynb` | Generate three-dimensional visualisations. |
| `sonde_flux_profiles.ipynb` | Plot vertical profiles of momentum flux. |
| `sonde_flux_wind_dir.ipynb` | Analyse momentum-flux direction relative to wind direction. |
| `sonde_height_profiles.ipynb` | Examine radiosonde altitude profiles. |
| `sonde_relative_humidity.ipynb` | Analyse relative-humidity profiles. |
| `sonde_trajectories_over_TOC.ipynb` | Visualise radiosonde trajectories. |

Notebooks may contain exploratory code, repeated calculations, and environment-specific file paths.

They should not currently be assumed to execute independently or in a particular order.

## 6. Dependencies

The python dependencies are outlined in the environment file: 

```text
teamx_environemnt.yml
```

The code uses scientific Python libraries for numerical calculations, interpolation, atmospheric data processing, and visualisation.

The dependency list has not yet been verified against a clean Python environment.

## 7. Running the Code

The repository currently contains Python source files and exploratory Jupyter notebooks.

A single entry point for executing the entire analysis has not yet been established.

The expected procedure is:

1. Ensure the required input NetCDF files are available in `data/`.
2. Install the dependencies from `requirements.txt`.
3. Run the relevant momentum-flux analysis code in `src/`.
4. Use the notebooks in `notebooks/` to inspect, compare, and visualise the resulting arrays.

Some file paths may need to be adjusted for the local environment.

The precise execution commands and required arguments will be documented once the source files have been reviewed.

## 8. Outputs (`results/`)

Selected figures are stored in:

```text
results/output_plots/
```

Current outputs include:

| File | Description |
|---|---|
| `Flux_mag_wind_dir.png` | Momentum-flux magnitude and wind-direction visualisation. |
| `Flux_magnitude_and_direction_iop_17_dpi1000.png` | Momentum-flux magnitude and direction for IOP 17. |
| `max_height_sondes.png` | Radiosonde maximum-height visualisation. |
| `trajectory_wind_component.png` | Radiosonde trajectory and wind-component visualisation. |

Additional numerical outputs are stored in the corresponding directories under `data/`.

## 9. Current Development Status

This repository contains research code that is being organised for sharing.

Documentation is incomplete in the following areas:

- Exact function inputs, outputs, and units.
- Dependencies between Python modules.
- Execution order and entry points.
- File-path configuration.
- Structure of the saved `.npz` archives.
- Reproducibility of individual notebooks.

These details will be added as the code is reviewed.

---

**Note:** This README is a preliminary technical description of the repository. File descriptions and execution relationships should be checked against the implementation.
