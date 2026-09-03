#!/usr/bin/env python3
"""Pure-xarray Vialard–Delecluse mixed-layer heat budget in W m-2."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import xarray as xr

RHO0 = 1035.0
CP = 3989.24495
EARTH_RADIUS = 6_371_000.0

def open_var(root: Path, stem: str, nmonths: int | None, chunks: int):
    path = next(root.glob(f"{stem}_Omon_*.nc"))
    ds = xr.open_dataset(path)
    chunk_map = {d: (-1 if d != "time" else chunks) for d in ds.dims}
    ds = ds.chunk(chunk_map)
    da = ds[stem]
    return da.isel(time=slice(0, nmonths)) if nmonths else da

def last_where(field, condition, zdim):
    """Deepest valid value satisfying condition, without eager 4-D indexing."""
    out = xr.full_like(field.isel({zdim: 0}, drop=True), np.nan)
    for k in range(field.sizes[zdim]):
        out = xr.where(condition.isel({zdim: k}), field.isel({zdim: k}), out)
    return out

def first_where(field, condition, zdim):
    """Shallowest valid value satisfying condition."""
    out = xr.full_like(field.isel({zdim: 0}, drop=True), np.nan)
    for k in range(field.sizes[zdim] - 1, -1, -1):
        out = xr.where(condition.isel({zdim: k}), field.isel({zdim: k}), out)
    return out

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, default=Path("../data/cm2"))
    p.add_argument("--output", default="ACCESS_CM2_HBA_python_Wm2.nc")
    p.add_argument("--nmonths", type=int, default=None, help="e.g. 600 for first 50 years")
    p.add_argument("--chunk-time", type=int, default=24)
    p.add_argument("--kv", type=float, default=1e-5, help="constant E2 diffusivity, m2 s-1")
    args = p.parse_args()

    get = lambda name: open_var(args.data_dir, name, args.nmonths, args.chunk_time)
    temp, mld, u, v, wt = get("temp"), get("mld"), get("u"), get("v"), get("wt")
    sw, coupler, pme, frazil = get("sw_heat"), get("sfc_hflux_coupler"), get("sfc_hflux_pme"), get("frazil_3d")
    zdim = "st_ocean"
    z = temp[zdim]

    # Infer cell interfaces from tracer-cell centres. ACCESS upper levels are
    # 5,15,... m, giving 0,10,... m interfaces.
    edges = np.empty(z.size + 1)
    zv = z.values.astype(float)
    edges[0] = max(0.0, zv[0] - 0.5 * (zv[1] - zv[0]))
    edges[1:-1] = 0.5 * (zv[:-1] + zv[1:])
    edges[-1] = zv[-1] + 0.5 * (zv[-1] - zv[-2])
    dz = xr.DataArray(np.diff(edges), dims=zdim, coords={zdim: z})
    bottom = xr.DataArray(edges[1:], dims=zdim, coords={zdim: z})
    inside = z < mld
    thick = dz.where(inside)
    h = thick.sum(zdim, skipna=True).where(inside.any(zdim)).rename("H_EXTN")
    tm = ((temp * thick).sum(zdim, skipna=True) / h).rename("T_MIX")

    # Real-time derivatives in SI seconds.
    sec = ((temp.time - temp.time[0]) / np.timedelta64(1, "s")).astype("float64")
    tm_sec = tm.assign_coords(time=sec)
    h_sec = h.assign_coords(time=sec)
    dtmp_dt = tm_sec.differentiate("time").assign_coords(time=temp.time)
    dh_dt = h_sec.differentiate("time").assign_coords(time=temp.time)
    tend = (RHO0 * CP * h * dtmp_dt).rename("TEND")

    # Surface forcing retained by the mixed layer.
    sw_below = sw.sum(zdim, skipna=True) - sw.where(inside).sum(zdim, skipna=True)
    shf = (coupler + pme + frazil.isel({zdim: 0}) - sw_below).rename("SHF")

    # Horizontal gradients on a regular lon/lat analysis grid. Differentiate
    # per degree, then convert degrees to metres locally.
    deg_to_rad = np.pi / 180.0
    dx_per_degree = EARTH_RADIUS * deg_to_rad * np.cos(np.deg2rad(temp.lat))
    dy_per_degree = EARTH_RADIUS * deg_to_rad
    dtdx = temp.differentiate("lon") / dx_per_degree
    dtdy = temp.differentiate("lat") / dy_per_degree
    ax = (-RHO0 * CP * (u * dtdx * thick).sum(zdim, skipna=True)).rename("AX")
    ay = (-RHO0 * CP * (v * dtdy * thick).sum(zdim, skipna=True)).rename("AY")

    # Temperature immediately below MLD and vertical velocity at its base.
    tb = first_where(temp, z >= mld, zdim)
    wdim = "sw_ocean" if "sw_ocean" in wt.dims else zdim
    wz = wt[wdim]
    wb = last_where(wt, wz <= mld, wdim)
    delta = tb - tm
    az = (RHO0 * CP * wb * delta).rename("AZ")
    e1 = (RHO0 * CP * xr.where(dh_dt > 0, dh_dt, 0) * delta).rename("E1")

    # Constant-Kv approximation at the first tracer point below MLD.
    dtdz = temp.differentiate(zdim)
    grad_b = first_where(dtdz, z >= mld, zdim)
    e2 = (RHO0 * CP * args.kv * grad_b).rename("E2")

    rhs = (shf + ax + ay + az + e1 + e2).rename("RHS")
    residual = (tend - rhs).rename("RESIDUAL")
    out = xr.Dataset({"T_MIX": tm, "H_EXTN": h, "TEND": tend, "SHF": shf,
                      "AX": ax, "AY": ay, "AZ": az, "E1": e1, "E2": e2,
                      "RHS": rhs, "RESIDUAL": residual})
    out.T_MIX.attrs.update(units="degC", long_name="mixed-layer mean temperature")
    out.H_EXTN.attrs.update(units="m", long_name="discrete mixed-layer thickness")
    for name in ("TEND", "SHF", "AX", "AY", "AZ", "E1", "E2", "RHS", "RESIDUAL"):
        out[name].attrs["units"] = "W m-2"
    out.attrs.update(
        equation="RESIDUAL = TEND - (SHF + AX + AY + AZ + E1 + E2)",
        rho0=RHO0, heat_capacity=CP, constant_Kv=args.kv,
        warning="E1 and constant-Kv E2 are reconstructed candidates, not verified one-to-one online diagnostics",
    )
    encoding = {name: {"zlib": True, "complevel": 1} for name in out.data_vars}
    out.to_netcdf(args.output, encoding=encoding)
    print(f"wrote {args.output}")

if __name__ == "__main__":
    main()
