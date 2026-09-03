"""Reconstruct mixed-layer heat storage and surface forcing in W m-2."""
from pathlib import Path

import numpy as np
import xarray as xr

RHO0 = 1035.0
CP = 3989.24495


def _open(data_dir, stem):
    path = Path(data_dir) / f"{stem}_Omon_bj594_piControl_r1i1p1_pac.nc"
    return xr.open_dataset(path)


def calculate_mixed_layer_budget(data_dir="../data/cm2",
                                 output_filename="mixed_layer_budget_output.nc"):
    temp = _open(data_dir, "temp")["temp"]
    mld = _open(data_dir, "mld")["mld"]
    sw = _open(data_dir, "sw_heat")["sw_heat"]
    coupler = _open(data_dir, "sfc_hflux_coupler")["sfc_hflux_coupler"]
    pme = _open(data_dir, "sfc_hflux_pme")["sfc_hflux_pme"]
    frazil = _open(data_dir, "frazil_3d")["frazil_3d"].isel(st_ocean=0)
    vdiff = _open(data_dir, "temp_vdiffuse_sbc")["temp_vdiffuse_sbc"]

    zdim = "st_ocean"
    z = temp[zdim]
    # Match the six-level Ferret mask: centres are used as nominal cell tops.
    cell_top = z
    cell_bot = z.shift({zdim: -1})
    last_bottom = z[-1] + (z[-1] - z[-2])
    cell_bot = cell_bot.fillna(last_bottom)
    inside = cell_top < mld
    h_extn = cell_bot.where(inside).max(zdim, skipna=True)

    # Ferret Z=@AVE is thickness weighted; all retained input levels are 10 m.
    dz = cell_bot - cell_top
    t_mix = (temp.where(inside) * dz).sum(zdim, skipna=True) / dz.where(inside).sum(zdim, skipna=True)

    # Differentiate using the real decoded time coordinate, expressed in seconds.
    seconds = ((temp.time - temp.time[0]) / np.timedelta64(1, "s")).astype("float64")
    t_seconds = t_mix.assign_coords(time=seconds)
    dtdt = t_seconds.differentiate("time").assign_coords(time=temp.time)
    tend = RHO0 * CP * h_extn * dtdt

    # Same Qnet definition as the Ferret reference calculation.
    sw_below_mld = sw.sum(zdim, skipna=True) - sw.where(inside).sum(zdim, skipna=True)
    shf_qnet = coupler + pme + frazil - sw_below_mld

    # Independent model diagnostic. It is non-zero mainly at the surface and is
    # retained separately; it must not be silently substituted for Qnet.
    shf_vdiff = vdiff.where(inside).sum(zdim, skipna=True)

    out = xr.Dataset({
        "T_MIX": t_mix,
        "H_EXTN": h_extn,
        "TEND": tend,
        "SHF": shf_qnet,
        "SHF_VDIFF": shf_vdiff,
        "SHF_DIFFERENCE": shf_qnet - shf_vdiff,
    })
    out["T_MIX"].attrs.update(units="degC", long_name="mixed-layer mean temperature")
    out["H_EXTN"].attrs.update(units="m", long_name="discrete extended mixed-layer depth")
    for name in ("TEND", "SHF", "SHF_VDIFF", "SHF_DIFFERENCE"):
        out[name].attrs["units"] = "W m-2"
    out["TEND"].attrs["long_name"] = "mixed-layer heat storage rho0 Cp H dT/dt"
    out["SHF"].attrs["long_name"] = "Qnet: coupler + PME + frazil - shortwave below MLD"
    out.attrs.update(rho0=RHO0, heat_capacity=CP,
                     warning="Reference SEF did not match Qnet or temp_vdiffuse_sbc during audit")
    out.to_netcdf(output_filename)
    print(f"Wrote {output_filename}")


if __name__ == "__main__":
    calculate_mixed_layer_budget()
