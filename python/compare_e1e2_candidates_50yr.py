#!/usr/bin/env python3
"""Test offline E1+E2 against model heat-budget remainder and process sums."""
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 600
RHO0, CP = 1035.0, 3989.24495
DATA = Path("../data/cm2")

def opened(stem):
    p = next(DATA.glob(f"{stem}_*.nc"))
    return xr.open_dataset(p).isel(time=slice(0, N))[stem]

b = xr.open_dataset("mixed_layer_budget_output.nc").isel(time=slice(0, N))
a = xr.open_dataset("ACCESS_CM2_adv_only_wm2_v8.nc").rename(
    {"TIME": "time", "LAT": "lat", "LON": "lon"}
).isel(time=slice(0, N))
adv = xr.open_dataset("temp_advection_Omon_bj594_piControl_pac_top50.nc").temp_advection.isel(time=slice(0, N))

time = xr.DataArray(b.time.values, dims="time")
lat = xr.DataArray(b.lat.values, dims="lat")
lon = xr.DataArray(b.lon.values, dims="lon")
a = a.assign_coords(time=time, lat=lat, lon=lon)
adv = adv.assign_coords(time=time, lat=lat, lon=lon)

def align_time(x):
    return x.assign_coords(time=time)

def mld_sum(x):
    x = align_time(x)
    if "st_ocean" in x.dims:
        return x.where(x.st_ocean < b.H_EXTN).sum("st_ocean", skipna=True)
    return x

# JNL formulas are temperature rates; rho Cp H converts them to W m-2.
offline_e12 = (RHO0 * CP * b.H_EXTN * (a.E_T1 + a.E_T2)).rename("OFFLINE_E1_PLUS_E2")

model_tendency = mld_sum(opened("temp_tendency"))
model_advection = adv.where(adv.st_ocean < b.H_EXTN).sum("st_ocean", skipna=True)
sw = align_time(opened("sw_heat"))
sw_below = sw.sum("st_ocean", skipna=True) - sw.where(sw.st_ocean < b.H_EXTN).sum("st_ocean", skipna=True)
model_sef = b.SHF_VDIFF - sw_below
model_remainder = (model_tendency - model_sef - model_advection).rename("MODEL_REMAINDER")

mixing_names = ["temp_vdiffuse_diff_cbt", "temp_nonlocal_KPP", "temp_vdiffuse_k33", "temp_sigma_diff"]
other_names = ["neutral_diffusion_temp", "neutral_gm_temp", "temp_submeso", "mixdownslope_temp", "temp_rivermix"]
mixing_sum = sum(mld_sum(opened(n)) for n in mixing_names).rename("MODEL_VERTICAL_MIXING_SUM")
other_sum = sum(mld_sum(opened(n)) for n in other_names)
eta = mld_sum(opened("temp_eta_smooth"))
explicit_sum = (mixing_sum + other_sum + eta).rename("MODEL_EXPLICIT_OTHER_SUM")

w = np.cos(np.deg2rad(b.lat))
def avg(x): return x.weighted(w).mean(("lat", "lon"), skipna=True)
def stat(ref, test):
    x, y = np.asarray(avg(ref)), np.asarray(avg(test)); ok=np.isfinite(x)&np.isfinite(y)
    return (float(np.corrcoef(x[ok],y[ok])[0,1]), float(np.mean(np.abs(y[ok]-x[ok]))),
            float(np.mean(y[ok]-x[ok])), float(np.sqrt(np.mean((y[ok]-x[ok])**2))))

refs = [("Tendency - SEF - advection", model_remainder),
        ("Explicit archived process sum", explicit_sum),
        ("Vertical-mixing candidates only", mixing_sum)]
series = {"OFFLINE_E1_PLUS_E2": avg(offline_e12)}
fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
for ax, (label, ref) in zip(axes, refs):
    r, mae, bias, rms = stat(ref, offline_e12)
    key = label.upper().replace(" ", "_").replace("-", "_")
    series[key] = avg(ref)
    print(f"{label}: r={r:.6g}, MAE={mae:.6g}, bias(offline-ref)={bias:.6g}, RMS={rms:.6g} W/m2")
    ax.plot(time, avg(ref), color="0.2", lw=1.1, label=label)
    ax.plot(time, avg(offline_e12), "--", color="crimson", lw=1.0, label="Offline E1+E2")
    ax.set_title(f"r={r:.3f}, MAE={mae:.2f}, bias={bias:.2f}, RMS={rms:.2f} W m$^{{-2}}$")
    ax.set_ylabel("W m$^{-2}$"); ax.grid(ls=":", alpha=.4); ax.legend(loc="upper right")
axes[-1].set_xlabel("Time")
fig.suptitle("E1+E2 attribution tests — first 50 years", weight="bold")
fig.tight_layout(); fig.savefig("ACCESS_CM2_E1_E2_candidates_50yr.png", dpi=200, bbox_inches="tight")
xr.Dataset(series).assign_attrs(
    note="Area-weighted time series; hypotheses, not one-to-one diagnostic identities"
).to_netcdf("ACCESS_CM2_E1_E2_candidates_50yr.nc")
print("wrote ACCESS_CM2_E1_E2_candidates_50yr.nc/.png")
