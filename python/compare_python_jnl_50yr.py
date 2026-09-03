#!/usr/bin/env python3
"""Compare the pure-Python and Ferret/JNL budgets for the first 50 years."""
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N=600; RHO0=1035.; CP=3989.24495; SEC_MONTH=30.*86400.
py=xr.open_dataset("ACCESS_CM2_HBA_python_50yr_Wm2.nc").isel(time=slice(0,N))
base=xr.open_dataset("ACCESS_CM2_ml_heatb_adv_wm2_v8.nc").rename(
    {"TIME":"time","LAT":"lat","LON":"lon"}).isel(time=slice(0,N))
adv=xr.open_dataset("ACCESS_CM2_adv_only_wm2_v8.nc").rename(
    {"TIME":"time","LAT":"lat","LON":"lon"}).isel(time=slice(0,N))
common_time=xr.DataArray(py.time.values,dims="time")
for ds in (base,adv): ds.coords["time"]=common_time
adv=adv.assign_coords(lat=py.lat,lon=py.lon)
base=base.assign_coords(lat=py.lat,lon=py.lon)
rate_factor=RHO0*CP*py.H_EXTN
jnl={
 "TEND":base.DT_DT*rate_factor,
 "SHF":base.SEF*rate_factor,
 "AX":adv.AX,"AY":adv.AY,"AZ":adv.AZ,
 "E1":adv.E_T1*rate_factor,"E2":adv.E_T2*rate_factor,
}
jnl["RHS"]=sum(jnl[n] for n in ("SHF","AX","AY","AZ","E1","E2"))
jnl["RESIDUAL"]=jnl["TEND"]-jnl["RHS"]
w=np.cos(np.deg2rad(py.lat))
def avg(x): return x.weighted(w).mean(("lat","lon"),skipna=True)
def stats(x,y):
 x,y=np.asarray(x),np.asarray(y); ok=np.isfinite(x)&np.isfinite(y)
 return np.mean(np.abs(y[ok]-x[ok])),np.mean(y[ok]-x[ok]),np.corrcoef(x[ok],y[ok])[0,1]
names=["TEND","SHF","AX","AY","AZ","E1","E2","RESIDUAL"]
out={}
fig,axes=plt.subplots(len(names),1,figsize=(14,22),sharex=True)
for ax,n in zip(axes,names):
 js,ps=avg(jnl[n]),avg(py[n]); mae,bias,r=stats(js,ps)
 print(f"{n}: MAE={mae:.6g}, bias(Python-JNL)={bias:.6g} W/m2, r={r:.6g}")
 out[f"{n}_JNL"]=js; out[f"{n}_PYTHON"]=ps
 ax.plot(js.time,js,color="0.2",lw=1.1,label="JNL")
 ax.plot(ps.time,ps,"--",color="crimson",lw=1,label="Pure Python")
 ax.axhline(0,color="0.7",lw=.6); ax.grid(ls=":",alpha=.4)
 ax.set_ylabel("W m$^{-2}$"); ax.set_title(f"{n}: MAE={mae:.2f}, bias={bias:.2f}, r={r:.3f}")
axes[0].legend(loc="upper right"); axes[-1].set_xlabel("Time")
fig.suptitle("Vialard–Delecluse budget: JNL vs pure Python — first 50 years",weight="bold")
fig.tight_layout(); fig.savefig("ACCESS_CM2_Python_vs_JNL_50yr_Wm2.png",dpi=200,bbox_inches="tight")
xr.Dataset(out).assign_attrs(note="Area-weighted comparison; all terms W m-2").to_netcdf("ACCESS_CM2_Python_vs_JNL_50yr_Wm2.nc")
print("wrote ACCESS_CM2_Python_vs_JNL_50yr_Wm2.nc/.png")
