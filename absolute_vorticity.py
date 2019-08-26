#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Sat Aug 24 14:33:16 EDT 2019
#import os, os.path, sys
import xarray as xr, numpy as np, pandas as pd
#import matplotlib.pyplot as plt
from numpy import pi, sin

def absolute_vorticity(vort850, lat):
    '''calculate absolute vorticity at 850hPa given relative vorticity (s**-1) and grid latitudes (degN)'''
    eta = vort850 + 2 * (2*pi/24/3600) * sin(lat*pi/180)
    eta.attrs['long_name'] = '850hPa absolute vorticity'
    eta.attrs['units'] = 's**-1'

    return eta

if __name__ == '__main__':
    ifile = '/tigress/wenchang/MODEL_OUT/CTL1860_noleap_tigercpu_intelmpi_18_576PE/POSTP/10000101.atmos_month.nc'
    ds = xr.open_dataset(ifile).sel(pfull=slice(100, None))
    is_ocean = ds.land_mask.load() < 0.1
    vort850 = ds.vort850
    lat = ds.grid_yt
    eta = absolute_vorticity(vort850, lat)
