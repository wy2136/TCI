#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Fri Aug 23 11:53:11 EDT 2019
#import os, os.path, sys
#import xarray as xr, numpy as np, pandas as pd
#import matplotlib.pyplot as plt
import xarray as xr
from numpy import exp, log

# physical parameters
Lv = 2.555e6 # J/kg
g = 9.81 #m/s**-2
c_p = 1005.7 # J/K/kg
Rd = 287 # J/K/kg
Rv = 461 # J/K/kg
epsilon = Rd/Rv
T0 = 273.15 # K

def saturated_vapor_pressure(T):
    '''calculated saturated water vapor pressure (in Pa) given temperature T (in Kelvin)'''
    return 610.94*exp(17.625*(T-T0)/(T-T0+243.04))

def mixing_ratio(p, e):
    '''calculate mixing ration given air pressure (p) and water vapor pressure (e)'''
    return epsilon*e/(p-e)

def mixing_ratio_by_q(q):
    '''calculate mixing ratio given specific humidity q'''
    return q/(1-q)

def vapor_pressure_by_mixing_ratio(p, r_v):
    '''calculate water vapor pressure given air pressure and mixing ratio'''
    return p*r_v/(epsilon+r_v)

def moist_entropy(T, p, RH=None, q=None):
    '''calculate moist entropy given air temperature (T), pressure (p) and relative humidity (RH, 0-1) or specific humidity (q).
    The equation is: s = c_p*log(T) - Rd*log(p_d) + Lv*r_v/T - Rv*r_v*log(RH)'''
    if RH is None:
        assert q is not None, 'at least one of the two variables must be specified: relative humidity/specific humidity'
        r_v = mixing_ratio_by_q(q)
        e = vapor_pressure_by_mixing_ratio(p, r_v)
        RH = e/saturated_vapor_pressure(T)
    else: # RH is from input directly
        e = saturated_vapor_pressure(T) * RH
        r_v = mixing_ratio(p, e)
    s = c_p*log(T) - Rd*log(p-e) + Lv*r_v/T - Rv*r_v*log(RH)
    s.attrs['long_name'] = 'moist entropy'
    s.attrs['units'] = 'J/K/kg'

    return s

def entropy_deficit(sst, slp, Tb, RHb, p_m, Tm, RHm):
    '''calculate entropy deficity defined in Tang and Emanuel, 2012.
    sst: sea surface temperature (in Kelvin);
    slp: sea level pressure (in Pa);
    Tb: boundary layer air temperature;
    RHb: boundary layer relative humidity (0-1);
    p_m: middle troposphere pressure level (usually 6e4 Pa);
    Tm: middle troposphere air temperature;
    RHm: middle troposphere relative humidity (0-1).'''
    s_sst_star = moist_entropy(T=sst, p=slp, RH=1)
    s_b        = moist_entropy(T=Tb, p=slp, RH=RHb)
    s_m_star   = moist_entropy(T=Tm, p=p_m, RH=1)
    s_m        = moist_entropy(T=Tm, p=p_m, RH=RHm)

    chi = (s_m_star - s_m)/(s_sst_star - s_b).pipe(lambda x: x.where(x>0)) # exclude values <= 0 
    chi.attrs['long_name'] = 'entropy deficit'

    return chi

if __name__ == '__main__':
    ifile = '/tigress/wenchang/MODEL_OUT/CTL1860_noleap_tigercpu_intelmpi_18_576PE/POSTP/10000101.atmos_month.nc'
    ofile = 'test.nc'
    ds = xr.open_dataset(ifile)
    land_mask = ds.land_mask.load()
    is_ocean = land_mask < 0.1
    
    p_m = 6e4 # 600hPa
    chi = entropy_deficit(sst=ds.t_surf,
        slp=ds.ps,
        Tb=ds.t_ref,
        RHb=ds.rh_ref/100,
        p_m=p_m,
        Tm=ds.temp.interp(pfull=p_m/100).drop('pfull'),
        RHm=ds.rh.interp(pfull=p_m/100).drop('pfull')/100)
