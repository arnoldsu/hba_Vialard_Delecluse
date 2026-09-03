# Vialard–Delecluse mixed-layer heat budget

This directory contains a reproducible ACCESS-CM2/OM2 mixed-layer heat-budget workflow based on the Vialard and Delecluse decomposition. Calculated heat-budget terms and comparisons are expressed in **W m-2**.

## Current status

The first 50 years (600 monthly records) have been tested for the Pacific subset used in this project.

| Offline term | Model/reference term | MAE (W m-2) | Correlation |
|---|---|---:|---:|
| `TEND` | converted online `DT_DT` | 0.187 | 0.99994 |
| `AX + AY + AZ` | model `temp_advection`, top 5 levels | 3.357 | 0.772 |
| `SHF` | `temp_vdiffuse_sbc - SW_below_MLD` | 0.961 | 0.99947 |

The signs of tendency, surface forcing and total advection are supported by the model diagnostics. Individual AX, AY and AZ cannot be validated separately with `temp_advection`, which contains only their sum.

## Budget equation

The convention is positive heating of the mixed layer:

```text
TEND = SHF + AX + AY + AZ + E1 + E2 + RESIDUAL
```

with

```text
TEND = rho0 Cp H d(Tm)/dt
AX   = -rho0 Cp integral[u dT/dx dz]
AY   = -rho0 Cp integral[v dT/dy dz]
AZ   =  rho0 Cp wb (Tb - Tm)
E1   =  rho0 Cp max(dH/dt, 0) (Tb - Tm)
E2   =  rho0 Cp Kv (dT/dz at H)
```

`Tm` is mixed-layer mean temperature, `Tb` is the first temperature below the mixed layer, `H` is the discrete extended MLD and `wb` is vertical velocity at the mixed-layer base.

## Surface heat flux

The offline mixed-layer surface forcing is

```text
SHF = sfc_hflux_coupler + sfc_hflux_pme + frazil_3d(surface) - SW_below_MLD
```

and

```text
sfc_hflux_coupler = sens_heat + evap_heat + lw_heat + swflx
```

`temp_vdiffuse_sbc` represents the full surface-boundary contribution. Therefore the like-for-like SEF reference is `temp_vdiffuse_sbc - SW_below_MLD`. The earlier approximately 9 W m-2 offset was principally penetrative shortwave, not a missing turbulent or radiative surface component.

## E1 and E2 limitation

E1 is a moving-boundary entrainment term and normally has no single online ACCESS/OM2 diagnostic. It must be reconstructed from MLD and temperature.

### Working assumption for model-output attribution

E1 and E2 are mathematical terms in the Vialard–Delecluse decomposition; they do not necessarily have a one-to-one correspondence with a single ACCESS/OM2 diagnostic. Our working assumption is therefore:

```text
E1 may be represented by one model heat-budget term,
or by the sum/residual of processes that move and mix across the MLD.

E2 may be represented by one vertical-mixing diagnostic,
or by the sum of several vertical diffusion/KPP diagnostics.
```

Candidate attribution is:

| Offline term | Possible ACCESS/OM2 representation |
|---|---|
| E1: moving-MLD entrainment | Direct reconstruction from `mld` and `temp`; alternatively part of the residual after all fixed-volume online heat-budget terms are summed |
| E2: mixing through the MLD base | One of `temp_vdiffuse_diff_cbt`, `temp_nonlocal_KPP`, or `temp_vdiffuse_k33`; more plausibly their MLD-integrated sum, with `temp_sigma_diff` tested as an additional small contribution |
| E1 + E2 | The difference between model `temp_tendency` and the sum of surface forcing, advection and all other explicitly available heating terms, provided every term uses the same sign, MLD mask and W m-2 units |

Thus the following alternatives should be compared quantitatively rather than assumed equivalent:

```text
E2a = sum_MLD(temp_vdiffuse_diff_cbt)
E2b = sum_MLD(temp_nonlocal_KPP)
E2c = sum_MLD(temp_vdiffuse_k33)
E2sum = E2a + E2b + E2c [+ sum_MLD(temp_sigma_diff)]

E12_residual = temp_tendency
             - SHF
             - total_advection
             - other explicitly diagnosed heat-budget processes
```

If `E1 + E2` resembles `E12_residual` but neither E1 nor E2 resembles an individual online field, that supports interpretation as a combined collection of model parameterized processes rather than two uniquely archived variables. A good correlation alone is insufficient: mean, MAE, sign, seasonal cycle and closure residual must also improve.

The first-50-year tests give:

| E1+E2 comparison | Correlation | MAE (W m-2) | Offline-reference bias (W m-2) |
|---|---:|---:|---:|
| `temp_tendency - SEF - temp_advection` | 0.525 | 13.92 | +13.92 |
| Sum of all available archived non-surface/non-advection processes | 0.534 | 9.19 | +9.19 |
| Vertical-mixing candidates only | 0.493 | 12.12 | +12.12 |

These moderate correlations support investigating a combined-process interpretation, but the large biases show that none of these sums is presently equivalent to offline E1+E2. Likely causes include moving-versus-fixed control volumes, incomplete archived processes, the constant-Kv E2 approximation, remapping and MLD-boundary discretization.

The JNL E2 uses constant `Kv = 1e-5 m2 s-1`. This is only a simple approximation to ACCESS-CM2 vertical mixing. Candidate model diagnostics are:

```text
temp_vdiffuse_diff_cbt
temp_nonlocal_KPP
temp_vdiffuse_k33
temp_sigma_diff
```

An MLD integral of the first three is a useful online vertical-mixing comparison, but it is not mathematically identical to constant-Kv E2. In the present 50-year test its correlation with the offline E2 is about 0.54. E1 and E2 must therefore remain labelled as reconstructed/candidate terms rather than fully validated diagnostics.

Important unit point: JNL E1/E2 formulas generate temperature tendency rates. Their W m-2 form is obtained with `rho0 Cp H` after confirming the time unit used by Ferret. Do not blindly apply a second seconds-per-month division.

## Files

- `jnl/hba_vialard_delecluse_wm2.jnl`: cleaned Ferret/JNL implementation used for AX, AY, AZ and experimental E1/E2.
- `python/calculate_heat_budget.py`: complete pure-Python replacement for the JNL; writes TEND, SHF, AX, AY, AZ, E1, E2, RHS and residual in W m-2.
- `python/calculate_tendency_shf.py`: offline tendency and surface-forcing reconstruction.
- `python/remap_temp_advection_top50.py`: memory-bounded nearest-neighbour diagnostic remap of the top five native levels.
- `python/compare_offline_reference_50yr.py`: three-panel 50-year validation plot and NetCDF output.
- `python/compare_e1e2_candidates_50yr.py`: compares offline E1+E2 with the model remainder, all explicit archived processes and vertical-mixing candidates.
- `python/compare_python_jnl_50yr.py`: final first-50-year, seven-term plus residual comparison between pure Python and JNL.
- `python/prepare_data.py`: input manifest checker and optional symlink creator.
- `examples/run_50yr.pbs`: PBS example. Adjust the project/module names for the current Gadi environment.

No large model data are stored in this directory.

## Required input data

The data directory is not sufficient by itself: it must contain the following NetCDF files and variables. Filenames may have a different experiment/member/date suffix, but the prefix and internal variable name must match unless the scripts are adapted.

### Required for the complete offline calculation

| Filename pattern | Variable | Expected dimensions | Expected units | Used for |
|---|---|---|---|---|
| `temp_Omon_*.nc` | `temp` | `time, st_ocean, lat, lon` | degC | Tm, Tb, TEND, temperature gradients, E1 and E2 |
| `mld_Omon_*.nc` | `mld` | `time, lat, lon` | m | Time-dependent mixed-layer mask and H |
| `u_Omon_*.nc` | `u` | `time, st_ocean, lat, lon` | m s-1 | AX |
| `v_Omon_*.nc` | `v` | `time, st_ocean, lat, lon` | m s-1 | AY |
| `wt_Omon_*.nc` | `wt` | `time, sw_ocean, lat, lon` (or compatible vertical dimension) | m s-1 | AZ at the MLD base |
| `sw_heat_Omon_*.nc` | `sw_heat` | `time, st_ocean, lat, lon` | W m-2 per layer | Shortwave penetrating below MLD |
| `sfc_hflux_coupler_Omon_*.nc` | `sfc_hflux_coupler` | `time, lat, lon` | W m-2 | Coupled sensible + latent + LW + SW flux |
| `sfc_hflux_pme_Omon_*.nc` | `sfc_hflux_pme` | `time, lat, lon` | W m-2 | Heat carried by precipitation-minus-evaporation mass transfer |
| `frazil_3d_Omon_*.nc` | `frazil_3d` | `time, st_ocean, lat, lon` | W m-2 | Frazil/sea-ice heat contribution |

The pure-Python calculation therefore reads `temp`, `mld`, `u`, `v`, `wt`, `sw_heat`, `sfc_hflux_coupler`, `sfc_hflux_pme` and `frazil_3d`. It does not need `temp_advection` or JNL output to calculate the budget.

### Required only for validation/comparison

| Filename/variable | Purpose |
|---|---|
| `temp_vdiffuse_sbc_Omon_*.nc` / `temp_vdiffuse_sbc` | Reference for the surface-boundary heat contribution; subtract `SW_below_MLD` for a like-for-like SEF reference |
| `temp_tendency_Omon_*.nc` / `temp_tendency` | Model heat-tendency reference and E1+E2 remainder calculation |
| Native or prepared `temp_advection` | Reference for total `AX+AY+AZ`; the current prepared file is `temp_advection_Omon_bj594_piControl_pac_top50.nc` |
| JNL outputs `ACCESS_CM2_ml_heatb_adv_wm2_v8.nc` and `ACCESS_CM2_adv_only_wm2_v8.nc` | Only required by the Python-versus-JNL comparison script |

### Optional process diagnostics for E1/E2 attribution

| Filename pattern / variable | Process |
|---|---|
| `temp_vdiffuse_diff_cbt_Omon_*.nc` | Convective/background vertical diffusion |
| `temp_nonlocal_KPP_Omon_*.nc` | Non-local KPP heating |
| `temp_vdiffuse_k33_Omon_*.nc` | K33 vertical diffusion |
| `temp_sigma_diff_Omon_*.nc` | Sigma-coordinate diffusion |
| `neutral_diffusion_temp_Omon_*.nc` | Neutral diffusion |
| `neutral_gm_temp_Omon_*.nc` | GM stirring |
| `temp_submeso_Omon_*.nc` | Submesoscale tendency |
| `mixdownslope_temp_Omon_*.nc` | Downslope mixing |
| `temp_rivermix_Omon_*.nc` | River mixing |
| `temp_eta_smooth_Omon_*.nc` | Surface/eta smoothing |

All calculation inputs must describe the same monthly records and compatible grids. If u, v, w and tracer temperature remain on different native staggered grids, they must be placed consistently at tracer points or the native-grid flux-divergence method should be used. Do not bilinearly remap extensive face transports in watts before taking their divergence.

## Input preparation

Validate the existing ACCESS-CM2 directory:

```bash
python hba_Vialard_Delecluse/python/prepare_data.py \
  --source /g/data/p66/ars599/work_budget/data/cm2
```

The same variable names can be supplied from OM2 output. Coordinate names/grid staggering may differ, so the JNL dataset paths and any remapping must be adapted. Native-grid heat-transport divergence is preferred over differentiating already-remapped extensive face transports.

## Run with pure Python

Run commands from the parent project directory (`ba_26062026_jk`), because the defaults use `../data/cm2`:

```bash
# Check that all required input variables exist.
python hba_Vialard_Delecluse/python/prepare_data.py \
  --source /g/data/p66/ars599/work_budget/data/cm2

# Calculate the first 50 years (600 monthly records).
python hba_Vialard_Delecluse/python/calculate_heat_budget.py \
  --data-dir ../data/cm2 --nmonths 600 \
  --output ACCESS_CM2_HBA_python_50yr_Wm2.nc
```

For the complete record, omit `--nmonths`:

```bash
python hba_Vialard_Delecluse/python/calculate_heat_budget.py \
  --data-dir ../data/cm2 \
  --output ACCESS_CM2_HBA_python_full_Wm2.nc
```

This is the recommended JNL-free workflow. It reads `temp`, `mld`, `u`, `v`, `wt` and the surface-flux inputs directly; it does not use JNL AX/AY/AZ as calculation inputs.

## Run with Ferret/JNL

On Gadi, from the same parent project directory:

```bash
module purge
module load ferret
ferret -nojnl -memsize 4000 \
  -script hba_Vialard_Delecluse/jnl/hba_vialard_delecluse_wm2.jnl
```

The JNL currently reads `../data/cm2` and processes the complete time range found in `temp`. Its output name is set by the `ncfl` symbol near the beginning of the JNL. Change that symbol if an existing output must be preserved.

## Run through PBS

Check the account, storage and modules in the example, then submit it from the project directory:

```bash
qsub hba_Vialard_Delecluse/examples/run_50yr.pbs
qstat -u "$USER"
```

## Generate comparisons

After the required calculation and reference files exist:

```bash
# Offline Python terms against model/reference diagnostics.
python hba_Vialard_Delecluse/python/compare_offline_reference_50yr.py

# Offline E1+E2 attribution tests.
python hba_Vialard_Delecluse/python/compare_e1e2_candidates_50yr.py

# Final pure-Python versus JNL comparison.
python hba_Vialard_Delecluse/python/compare_python_jnl_50yr.py
```

E1 and constant-Kv E2 remain hypotheses even though both implementations calculate them completely.

### Pure Python versus JNL (first 50 years)

| Term | MAE (W m-2) | Correlation |
|---|---:|---:|
| TEND | 2.69 | 0.947 |
| SHF | 3.71 | 0.996 |
| AX | 0.44 | 0.999 |
| AY | 1.07 | 0.997 |
| AZ | 0.83 | 0.994 |
| E1 | 0.17 | 0.993 |
| E2 | 0.09 | 0.984 |
| residual | 2.25 | 0.847 |

These statistics use the actual JNL v8 output for TEND/SEF and the corrected JNL advection output for AX/AY/AZ/E1/E2. The close AX/AY/AZ agreement shows that the pure-Python metric derivatives and MLD integration reproduce the corrected JNL calculation well. E1/E2 agree numerically with the JNL definitions but remain physically unverified against a unique model diagnostic.

Expected comparison outputs are:

```text
ACCESS_CM2_offline_vs_reference_50yr_Wm2.nc
ACCESS_CM2_offline_vs_reference_50yr_Wm2.png
```

The advection reference currently uses a nearest-wet-cell diagnostic remap and native levels centred at 5, 15, 25, 35 and 45 m. It is suitable for sign, magnitude and temporal validation but is not a conservative remap.

## Recommended next step

For a final closed ACCESS/OM2 budget, compute all online diagnostic convergences on the native ocean grid first, integrate them through the same MLD mask, and only then conservatively remap the resulting two-dimensional W m-2 fields. This is especially important for advection and vertical mixing.
