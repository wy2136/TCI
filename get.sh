#!/usr/bin/env bash
# Wenchang Yang (wenchang@princeton.edu)
# Mon Aug 26 13:11:09 EDT 2019
idir=/tigress/wenchang/MODEL_OUT/CTL1860_noleap_tigercpu_intelmpi_18_576PE/POSTP
ncks -v t_surf,slp,land_mask,t_ref,rh_ref,temp,rh,sphum,ucomp,vcomp,vort850 $idir/10000101.atmos_month.nc 10000101.atmos_month.nc
