#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Sat Aug 24 15:02:58 EDT 2019
import os.path#, os, sys
import xarray as xr, numpy as np#, pandas as pd
#import matplotlib.pyplot as plt
from numpy import absolute, exp, log

from shared.entropy_deficit import entropy_deficit
from shared.potential_intensity import potential_intensity
from shared.wind_shear import wind_shear
from shared.absolute_vorticity import absolute_vorticity

def do_tci(ifile, odir=None):
    '''calculate TC indices (e.g. GPI, VI) and related variables given FLOR/HiRAM atmos_month output'''
    print('[input]:', ifile)
    if odir is None:
        odir = '.'
    #ifile = '/tigress/wenchang/MODEL_OUT/CTL1860_noleap_tigercpu_intelmpi_18_576PE/POSTP/10000101.atmos_month.nc'
    ibasename = os.path.basename(ifile)
    ds = xr.open_dataset(ifile)
    is_ocean = ds.land_mask.load() < 0.1

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
            sst=ds.t_surf,
            slp=ds.slp*100,
            Tb=ds.t_ref,
            RHb=ds.rh_ref/100,
            p_m=p_m,
            Tm=ds.temp.interp(pfull=p_m/100).drop('pfull'),
            RHm=ds.rh.interp(pfull=p_m/100).drop('pfull')/100
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
        truncate_and_reverse_plevels = lambda x: x.sel(pfull=slice(100, None)) \
            .isel(pfull=slice(-1, None, -1)) # 
        PI = potential_intensity(
            sst=ds.t_surf.where(is_ocean),
            slp=ds.slp.where(is_ocean)*100,
            p=ds.pfull.pipe(truncate_and_reverse_plevels),
            T=ds.temp.pipe(truncate_and_reverse_plevels).where(is_ocean),
            q=ds.sphum.pipe(truncate_and_reverse_plevels).where(is_ocean),
            dim_x='grid_xt', dim_y='grid_yt', dim_z='pfull'
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
            u850=ds.ucomp.interp(pfull=850),
            v850=ds.vcomp.interp(pfull=850),
            u200=ds.ucomp.interp(pfull=200),
            v200=ds.vcomp.interp(pfull=200)
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
            vort850=ds.vort850,
            lat=ds.grid_yt
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
        H = ds.rh.interp(pfull=600).drop('pfull')
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
    #ifile = '/tigress/wenchang/MODEL_OUT/CTL1860_noleap_tigercpu_intelmpi_18_576PE/POSTP/10000101.atmos_month.nc'
    ifile = '10000101.atmos_month.nc'
    do_tci(ifile)
