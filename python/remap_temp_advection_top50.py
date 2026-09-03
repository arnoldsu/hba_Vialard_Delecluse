"""Chunked nearest-wet-cell remap of native temp_advection, top 50 m only."""
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from scipy.spatial import cKDTree
SRC='/g/data/p66/ars599/work_budget/data/temp_advection_Omon_bj594_piControl_18500630-20141231.nc'
OUT='temp_advection_Omon_bj594_piControl_pac_top50.nc'
src=xr.open_dataset(SRC,decode_times=False); tgt=xr.open_dataset('../data/cm2/temp_Omon_bj594_piControl_r1i1p1_pac.nc',decode_times=False)
lon=(src.geolon_t.values+360)%360; lat=src.geolat_t.values
valid=np.isfinite(src.temp_advection.isel(time=0,st_ocean=0).values); iy,ix=np.where(valid)
def xyz(lo,la):
 lo=np.deg2rad(lo); la=np.deg2rad(la); return np.c_[np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)]
LON,LAT=np.meshgrid(tgt.lon.values,tgt.lat.values); tree=cKDTree(xyz(lon[valid],lat[valid])); _,q=tree.query(xyz(LON.ravel(),LAT.ravel()))
ny,nx=tgt.sizes['lat'],tgt.sizes['lon']; nt=src.sizes['time']; nz=5
ys,ye=int(iy[q].min()),int(iy[q].max())+1; xs,xe=int(ix[q].min()),int(ix[q].max())+1
print("native slab",ys,ye,xs,xe,flush=True)
with Dataset(OUT,'w') as nc:
 for n,size in [('time',nt),('st_ocean',nz),('lat',ny),('lon',nx)]: nc.createDimension(n,size)
 for n,data,dims in [('time',tgt.time.values,('time',)),('st_ocean',src.st_ocean.values[:nz],('st_ocean',)),('lat',tgt.lat.values,('lat',)),('lon',tgt.lon.values,('lon',))]: nc.createVariable(n,'f8',dims)[:]=data
 vout=nc.createVariable('temp_advection','f4',('time','st_ocean','lat','lon'),zlib=True,complevel=1,fill_value=np.nan); vout.units='W m-2'; vout.remapping='nearest native wet tracer cell'
 for start in range(0,nt,24):
  stop=min(start+24,nt)
  slab=src.temp_advection.isel(time=slice(start,stop),st_ocean=slice(0,nz),yt_ocean=slice(ys,ye),xt_ocean=slice(xs,xe)).values
  mapped=slab[:,:,iy[q]-ys,ix[q]-xs].reshape(stop-start,nz,ny,nx)
  vout[start:stop]=mapped
  print(start,stop,flush=True)
print('wrote',OUT)
