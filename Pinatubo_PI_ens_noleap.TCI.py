#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Tue Sep 17 11:11:47 EDT 2019
#import os.path, sys, os
#import xarray as xr, numpy as np, pandas as pd
#import matplotlib.pyplot as plt
import os, os.path, sys
from AMx import do_tci

expname = 'Pinatubo_PI_ens_noleap'
year_eruption = 1991
years = range(year_eruption,year_eruption+5)
ens = range(1,31)

exp_user = 'wenchang'
analysis_user = os.environ['USER']

for en in ens:
    print(f'en = {en:02d}' )
    idir = f'/tigress/{exp_user}/MODEL_OUT/{expname}/en{en:02d}'
    odir = f'/tigress/{analysis_user}/MODEL_OUT/{expname}/en{en:02d}/analysis_wy/TCI'
    if not os.path.exists(odir):
        os.makedirs(odir)
        print('[dir made]:', odir)
    for year in years:
        print(f'year = {year:04d}')
        ifile = os.path.join(idir, 'POSTP', f'{year:04d}0101.atmos_month.nc')
        do_tci(ifile, odir)
print('**done**')
