#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Sat Aug 24 15:02:58 EDT 2019
import os.path#, os, sys
import xarray as xr, numpy as np#, pandas as pd
#import matplotlib.pyplot as plt
from numpy import absolute, exp, log

from entropy_deficit import entropy_deficit
from potential_intensity import potential_intensity
from wind_shear import wind_shear
from absolute_vorticity import absolute_vorticity

def do_tci(year, odir=None):
    '''calculate TC indices (e.g. GPI, VI) and related variables given FLOR/HiRAM atmos_month output'''
    print('[year]:', year)
    if odir is None:
        odir = '.'
    ibasename = f'era5.monthly.{year}.nc' 

    # sst and ocean mask
    ifile = f'/tigress/wenchang/data/era5/raw/monthly/sea_surface_temperature/era5.sea_surface_temperature.monthly.{year:04d}.nc'
    sst = xr.open_dataarray(ifile)# units K
    is_ocean = sst.isel(time=0).drop('time').pipe(lambda x: x*0==0)

    # slp
    ifile = f'/tigress/wenchang/data/era5/raw/monthly/mean_sea_level_pressure/era5.mean_sea_level_pressure.monthly.{year:04d}.nc'
    slp = xr.open_dataarray(ifile) # units Pa
    
    # t2m
    ifile = f'/tigress/wenchang/data/era5/analysis/2m_temperature/monthly/era5.2m_temperature.monthly.{year:04d}.nc'
    t2m = xr.open_dataarray(ifile) # units K

    # RH2m
    ifile = f'/tigress/wenchang/data/era5/analysis/2m_relative_humidity/monthly/era5.2m_relative_humidity.monthly.{year:04d}.nc'
    RH2m = xr.open_dataarray(ifile) # units %
    
    # Ta
    ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/temperature/era5.temperature.monthly.{year:04d}.nc'
    Ta = xr.open_dataarray(ifile) # in K

    # RH
    ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/relative_humidity/era5.relative_humidity.monthly.{year:04d}.nc'
    RH = xr.open_dataarray(ifile) # in %

    # q
    ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/specific_humidity/era5.specific_humidity.monthly.{year:04d}.nc'
    q = xr.open_dataarray(ifile) # in kg/kg

    # u
    ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/u_component_of_wind/era5.u_component_of_wind.monthly.{year:04d}.nc'
    u = xr.open_dataarray(ifile) # in m/s

    # v 
    ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/v_component_of_wind/era5.v_component_of_wind.monthly.{year:04d}.nc'
    v = xr.open_dataarray(ifile) # in m/s

    # vorticity
    ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/vorticity/era5.vorticity.monthly.{year:04d}.nc'
    vort = xr.open_dataarray(ifile) # in s**-1


    # entropy deficit: (s_m_star - s_m)/(s_sst_star - s_b)
    print('entropy deficit')
    dname = 'chi'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        chi = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        p_m = 600*100 # Pa
        chi = entropy_deficit(
            sst=sst,
            slp=slp,
            Tb=t2m,
            RHb=RH2m/100,
            p_m=p_m,
            Tm=Ta.sel(level=p_m/100).drop('level'),
            RHm=RH.sel(level=p_m/100).drop('level')/100
            ).where(is_ocean)
        chi.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # potential intensity
    print('potential intensity')
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.PI.nc') )
    if os.path.exists(ofile):
        PI = xr.open_dataset(ofile)
        print('[opened]:', ofile)
    else:
        reverse_plevels = lambda x: x.isel(level=slice(-1, None, -1)) # 
        PI = potential_intensity(
            sst=sst,
            slp=slp.where(is_ocean),
            p=Ta.level.pipe(reverse_plevels),
            T=Ta.pipe(reverse_plevels).where(is_ocean),
            q=q.pipe(reverse_plevels).where(is_ocean),
            dim_x='longitude', dim_y='latitude', dim_z='level'
            )
        encoding = {dname:{'dtype': 'float32', 'zlib': True, 'complevel': 1} 
            for dname in ('pmin', 'vmax')}
        encoding['iflag'] = {'dtype': 'int32'}
        PI.to_netcdf(ofile, encoding=encoding, unlimited_dims='time')
        print('[saved]:', ofile)

    # wind shear: ( (u200-u850)**2 + (v200-v850)**2 )**0.5
    print('wind shear')
    dname = 'Vshear'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        Vshear = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        Vshear = wind_shear(
            u850=u.sel(level=850),
            v850=v.sel(level=850),
            u200=u.sel(level=200),
            v200=v.sel(level=200)
            )
        Vshear.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # ventilation index: Vshear * chi_m /V_PI
    print('ventilation index')
    dname = 'VI'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        VI = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        VI = Vshear*chi/PI.vmax.pipe(lambda x: x.where(x>0))
        VI.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # absolute vorticity at 850hPa
    print('absolute vorticity')
    dname = 'eta'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        eta = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        eta = absolute_vorticity(
            vort850=vort.sel(level=850),
            lat=vort.latitude
            )
        eta.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)


    # relative humidity at 600hPa in %
    print('relative humidity in %')
    dname = 'H'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        H = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        H = RH.sel(level=600).drop('level')
        H.attrs['long_name'] = '600hPa relative humidity'
        H.attrs['units'] = '%'
        H.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)
    
    # GPI (Emanuel and Nolan 2004): |10**5\eta|**(3/2) * (H/50)**3 * (Vpot/70)**3 * (1+0.1*Vshear)**(-2)
    print('GPI')
    dname = 'GPI'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        GPI = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        GPI = (1e5 * absolute(eta) )**(3/2) \
            * (H/50)**3 \
            * (PI.vmax/70)**3 \
            * (1+0.1*Vshear)**(-2)
        GPI.attrs['long_name'] = 'Genesis Potential Index'
        GPI.attrs['history'] = '|10**5\eta|**(3/2) * (H/50)**3 * (Vpot/70)**3 * (1+0.1*Vshear)**(-2)'
        GPI.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # GPI2010 (Emanuel 2010): |\eta|**3 * chi**(-4/3) * max((Vpot-35),0)**2 * (25+Vshear)**(-4)
    print('GPI2010')
    dname = 'GPI2010'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        GPI2010 = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        GPI2010 = absolute(eta)**3 \
            * chi.where(chi>0)**(-4/3) \
            * (PI.vmax - 35).clip(min=0)**2 \
            * (25 + Vshear)**(-4)
        GPI2010.attrs['long_name'] = 'Genesis Potential Index of Emanuel2010'
        GPI2010.attrs['history'] = '|\eta|**3 * chi**(-4/3) * max((Vpot-35),0)**2 * (25+Vshear)**(-4)'
        GPI2010.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)
    
if __name__ == '__main__':
    #year = 1979
    #do_tci(year, odir='/tigress/wenchang/data/era5/analysis/TCI/')
    years = range(1979, 2019)
    for year in years:
        do_tci(year, odir='/tigress/wenchang/data/era5/analysis/TCI/')
