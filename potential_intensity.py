#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Fri Aug 23 12:52:55 EDT 2019
#import os, os.path, sys
#import xarray as xr, numpy as np, pandas as pd
#import matplotlib.pyplot as plt
import numpy as np, xarray as xr

from pcmin import pcmin3

def potential_intensity(sst, slp, p, T, q, dim_x, dim_y, dim_z):
    '''xarray-wrapper of the FORTRAN module pcmin3_kflag.
    sst: sea surface temperature;
    slp: seal level pressure;
    p: pressure levels;
    T: temperature;
    q: specific humidity;
    xname: dim name along the x/lon direction;
    yname: dim name along the y/lat direction;
    zname: dim name along the z/p direction.
    '''
    r_v = q/(1-q) # specific humidity to mixing ratio
    pmin, vmax, iflag = xr.apply_ufunc(pcmin3,
        sst, slp, p, T, r_v,
        input_core_dims=[[dim_y, dim_x], [dim_y, dim_x], [dim_z], [dim_z, dim_y, dim_x], [dim_z, dim_y, dim_x]],
        output_core_dims=[[dim_y, dim_x], [dim_y, dim_x], []],
        vectorize=True)
    pmin.attrs['long_name'] = 'mininum central pressure'
    pmin.attrs['units'] = 'hPa'
    vmax.attrs['long_name'] = 'maximum surface wind speed'
    vmax.attrs['units'] = 'm/s'
    iflag = iflag.astype('int32')
    iflag.attrs['long_name'] = '1: OK; 0: no convergence; 2: CAPE routine failed.'
    PI = xr.Dataset(dict(pmin=pmin, vmax=vmax, iflag=iflag))
    
    return PI

if __name__ == '__main__':
    ifile = '/tigress/wenchang/MODEL_OUT/CTL1860_noleap_tigercpu_intelmpi_18_576PE/POSTP/10000101.atmos_month.nc'
    ds = xr.open_dataset(ifile).sel(pfull=slice(100, None))
    is_ocean = ds.land_mask.load() < 0.1
    sst = ds.t_surf.where(is_ocean)
    slp = ds.slp.where(is_ocean)*100
    i_reversed = slice(-1, None, -1)
    p = ds.pfull.isel(pfull=i_reversed)
    T = ds.temp.where(is_ocean).isel(pfull=i_reversed).load()
    q = ds.sphum.where(is_ocean).isel(pfull=i_reversed).load()
    PI = potential_intensity(sst, slp, p, T, q, dim_x='grid_xt', dim_y='grid_yt', dim_z='pfull')
    PI.to_netcdf('PI_test.nc', encoding={dname:{'zlib': True, 'complevel': 1} for dname in ('pmin', 'vmax')})
