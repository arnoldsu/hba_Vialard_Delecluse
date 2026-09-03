"""Compare the offline heat budget with model/reference diagnostics (first 50 years)."""
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NMONTH = 600
RHO0, CP, SEC_MONTH = 1035.0, 3989.24495, 30.0 * 86400.0

# Offline surface/tendency reconstruction and model surface-diffusion reference.
b = xr.open_dataset("mixed_layer_budget_output.nc").isel(time=slice(0, NMONTH))
# Original online mixed-layer tendency diagnostic (degC/month).
t = xr.open_dataset("ACCESS_CM2_ml_heatb_const_vdiff.nc").isel(time=slice(0, NMONTH))
# Corrected offline advection terms, already W m-2.
a = xr.open_dataset("ACCESS_CM2_adv_only_wm2_v8.nc").rename(
    {"TIME": "time", "LAT": "lat", "LON": "lon"}
).isel(time=slice(0, NMONTH))
# Native-model temp_advection remapped to the analysis grid; first five levels (5--45 m).
m = xr.open_dataset("temp_advection_Omon_bj594_piControl_pac_top50.nc").isel(
    time=slice(0, NMONTH)
)
sw = xr.open_dataset("../data/cm2/sw_heat_Omon_bj594_piControl_r1i1p1_pac.nc").sw_heat.isel(
    time=slice(0, NMONTH)
)

# These files contain the same consecutive monthly records but encode time differently.
common_time = xr.DataArray(b.time.values, dims="time", attrs=b.time.attrs)
for ds in (t, a, m):
    ds.coords["time"] = common_time
common_lat = xr.DataArray(b.lat.values, dims="lat", attrs=b.lat.attrs)
common_lon = xr.DataArray(b.lon.values, dims="lon", attrs=b.lon.attrs)
a = a.assign_coords(lat=common_lat, lon=common_lon)
m = m.assign_coords(lat=common_lat, lon=common_lon)

factor = RHO0 * CP * b.H_EXTN / SEC_MONTH
tendency_ref = (t.DT_DT * factor).rename("TENDENCY_REFERENCE")
tendency_off = b.TEND.rename("TENDENCY_OFFLINE")

# Each model temp_advection layer is W m-2. Select levels inside the MLD, capped by
# the available top-50-m reference, then vertically sum.
advection_ref = m.temp_advection.where(m.st_ocean < b.H_EXTN).sum(
    "st_ocean", skipna=True
).rename("ADVECTION_REFERENCE")
advection_off = (a.AX + a.AY + a.AZ).rename("ADVECTION_OFFLINE")

# _sbc contains the full surface boundary flux. Subtract shortwave that penetrates
# below the MLD so its definition matches the offline mixed-layer surface heating.
inside_sw = sw.st_ocean < b.H_EXTN
sw_below = sw.sum("st_ocean", skipna=True) - sw.where(inside_sw).sum("st_ocean", skipna=True)
shf_ref = (b.SHF_VDIFF - sw_below).rename("SHF_REFERENCE_SEF")
shf_off = b.SHF.rename("SHF_OFFLINE_SEF")

fields = [
    ("Tendency: DT_DT vs offline TEND", tendency_ref, tendency_off),
    ("Advection: model temp_advection vs AX+AY+AZ", advection_ref, advection_off),
    ("Surface heat flux: SEF (_sbc - SW below MLD) vs offline SHF", shf_ref, shf_off),
]
weights = np.cos(np.deg2rad(b.lat))

def spatial_mean(x):
    return x.weighted(weights).mean(("lat", "lon"), skipna=True)

def statistics(ref, off):
    x, y = np.asarray(ref), np.asarray(off)
    ok = np.isfinite(x) & np.isfinite(y)
    return {
        "mae": float(np.mean(np.abs(y[ok] - x[ok]))),
        "bias": float(np.mean(y[ok] - x[ok])),
        "corr": float(np.corrcoef(x[ok], y[ok])[0, 1]),
        "rms": float(np.sqrt(np.mean((y[ok] - x[ok]) ** 2))),
    }

out = xr.Dataset({
    tendency_ref.name: tendency_ref,
    tendency_off.name: tendency_off,
    advection_ref.name: advection_ref,
    advection_off.name: advection_off,
    shf_ref.name: shf_ref,
    shf_off.name: shf_off,
})
for var in out.data_vars:
    out[var].attrs["units"] = "W m-2"
out.attrs["period"] = "first 600 months (50 years)"
out.attrs["advection_note"] = "reference contains native model levels centred at 5,15,25,35,45 m"
out.to_netcdf(
    "ACCESS_CM2_offline_vs_reference_50yr_Wm2.nc",
    encoding={v: {"zlib": True, "complevel": 1} for v in out.data_vars},
)

fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
for ax, (title, ref, off) in zip(axes, fields):
    rs, os = spatial_mean(ref), spatial_mean(off)
    s = statistics(rs, os)
    print(f"{title}: MAE={s['mae']:.6g}, bias={s['bias']:.6g}, "
          f"RMS={s['rms']:.6g} W/m2, r={s['corr']:.6g}")
    ax.plot(rs.time, rs, color="0.2", lw=1.15, label="Model/reference")
    ax.plot(os.time, os, "--", color="crimson", lw=1.0, label="Offline")
    ax.axhline(0, color="0.7", lw=.7)
    ax.set_title(f"{title}\nMAE={s['mae']:.2f}, bias(off-ref)={s['bias']:.2f} W m$^{{-2}}$, r={s['corr']:.3f}")
    ax.set_ylabel("W m$^{-2}$")
    ax.grid(ls=":", alpha=.4)
axes[0].legend(loc="upper right")
axes[-1].set_xlabel("Time")
fig.suptitle("ACCESS-CM2 offline heat budget vs model diagnostics — first 50 years", weight="bold")
fig.tight_layout()
fig.savefig("ACCESS_CM2_offline_vs_reference_50yr_Wm2.png", dpi=200, bbox_inches="tight")
print("wrote ACCESS_CM2_offline_vs_reference_50yr_Wm2.nc/.png")
