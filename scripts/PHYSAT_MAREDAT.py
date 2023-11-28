#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Stephanie Anderson, Massachusetts Institute of Technology
Email: siander@mit.edu
    
This script examines a warming only climate change scenario 
(nutrients, circulation etc. are held constant)

    INPUT: 
        PHYSAT data (physat_monthly_mapped_degree.nc)
        MAREDAT data downloaded from http://www.pangaea.de/search?&q=maredat 
    
    OUTPUT: 
       Figure S7. PHYSAT
       Figure S8. MAREDAT
     
%=========================================================================
"""

#%% PHYSAT
import os 
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
from copy import copy
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

##### PHYSAT PFTs 
# Nano: 1
# Prochl: 2
# SLC: 3 (synechecoccus)
# Diatoms: 4
# Phaeocystis_like: 5
# Coccolithophorid(underestimate): 6

os.chdir("/Users/Stephanie/Desktop/MIT/Q10_Variability/Code & Datasets/")
pfts = (["Nanoeukaryotes", "Cyanobacteria", "Diatoms","Phaeocystis","Coccolithophores"])
colors1=["turquoise", "#ec3a25", "#026cb1", "lightgrey","orange"]
colors2=["white","turquoise", "#ec3a25", "#ec3a25", "#026cb1", "lightgrey", "orange"]
nodes = [0.0, 0.17, 0.34, 0.51,0.7,0.85,1.0]
cmap1 = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes,colors2)))
cmap1.set_bad(color = 'white', alpha = 1)

physat = 'data/physat_monthly_mapped_degree.nc'
physat = nc.Dataset(physat, 'r')

phyto = physat.variables['phyto']
lat = physat.variables['LATITUDE']
lon = physat.variables['LONGITUDE']

gridspec = dict(hspace= 0.001,bottom=0, height_ratios=[1, 1])
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6, 6),sharex='col', gridspec_kw=gridspec)

ax1 = plt.subplot(2,1,1,projection=ccrs.PlateCarree(central_longitude=0, globe=None))
im1 = ax1.imshow(phyto[0], extent=(-179,181,-90, 90), cmap=cmap1, vmin=0, vmax=6, 
                 origin='lower',  transform=ccrs.PlateCarree())
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax1.set_yticks([-50,0,50])
ax1.set_title(r'$\bf{(a)}$'+ ' January', loc='left', fontsize=14)

ax2 = plt.subplot(2,1,2,projection=ccrs.PlateCarree(central_longitude=0, globe=None))
im2 = ax2.imshow(phyto[6], extent=(-179,181,-90, 90), cmap=cmap1, vmin=0, vmax=6, 
                 origin='lower',  transform=ccrs.PlateCarree())
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax2.set_yticks([-50,0,50])
ax2.set_xticks([-100,0,100])
ax2.set_xlabel('Longitude', fontsize=12)
ax2.set_title(r'$\bf{(b)}$'+' July', loc='left', fontsize=14)

fig.text(0.04,0.44, 'Latitude', va='center', rotation='vertical', fontsize=12)

values = range(5)
patches = [mpatches.Patch(color=colors1[i], label=pfts[i]) for i in range(len(values))]
plt.legend(handles=patches, bbox_to_anchor=(1.25, 0.85), loc="lower center", borderaxespad=0., ncol=1, frameon=False)
plt.tight_layout()
#plt.savefig('figures/FigureS7.pdf', bbox_inches='tight', transparent=True)

#%% 
import numpy as np
# This data must be downloaded separately from http://www.pangaea.de/search?&q=maredat
maredat1 = '/Volumes/SABackup/MIT_q10/MAREDAT/MarEDat20130523Coccolithophores.nc' # https://doi.org/10.1594/PANGAEA.785092
maredat2 = '/Volumes/SABackup/MIT_q10/MAREDAT/MarEDat20120716Diatoms.nc' # https://doi.org/10.1594/PANGAEA.777384
maredat3 = '/Volumes/SABackup/MIT_q10/MAREDAT/MarEDat20130403Diazotrophs.nc' # https://doi.org/10.1594/PANGAEA.818214
maredat4 = '/Volumes/SABackup/MIT_q10/MAREDAT/MarEDat20111206Picophytoplankton.nc' # https://doi.org/10.1594/PANGAEA.777385
maredat1 = nc.Dataset(maredat1, 'r')
maredat2 = nc.Dataset(maredat2, 'r')
maredat3 = nc.Dataset(maredat3, 'r')
maredat4 = nc.Dataset(maredat4, 'r')

obs1 = maredat1.variables['NON_ZERO_BIOM']
obs2 = maredat2.variables['NON_ZERO_BIOM']
obs3 = maredat3.variables['NON_ZERO_BIOM']
obs4 = maredat4.variables['NON_ZERO_BIOM']
depth = maredat1.variables['DEPTH'][:]

# Calculate depth of each section
depth2 = np.empty([10])
depth_section = np.empty([10])
for x in range(10):
    depth2[x] = depth[x+1]
    depth_section[x] = depth2[x] - depth[x]

obs_yr1 = np.nansum(obs1, axis=0)
obs_int1 = np.empty([10,180,360])
obs_yr2 = np.nansum(obs2, axis=0)
obs_int2 = np.empty([10,180,360])
obs_yr3 = np.nansum(obs3, axis=0)
obs_int3 = np.empty([10,180,360])
obs_yr4 = np.nansum(obs4, axis=0)
obs_int4 = np.empty([10,180,360])

for x in range(10):
    obs_int1[x] = depth_section[x] * obs_yr1[x]*1000 # convert from L to m3
    obs_int2[x] = depth_section[x] * obs_yr2[x]*1000
    obs_int3[x] = depth_section[x] * obs_yr3[x]*1000
    obs_int4[x] = depth_section[x] * obs_yr4[x]*1000

cocco = np.nansum(obs_int1, axis=0)
diatoms = np.nansum(obs_int2, axis=0)
diazo = np.nansum(obs_int3, axis=0)
pico = np.nansum(obs_int4, axis=0)


#%% Figure S8

# Decreasing figure resolution to increase point sizes
cocco2 = cocco.reshape(45,4,90,4).sum(axis=1).sum(axis=2)
diatoms2 = diatoms.reshape(45,4,90,4).sum(axis=1).sum(axis=2)
pico2 = pico.reshape(45,4,90,4).sum(axis=1).sum(axis=2)
diazo2 = diazo.reshape(45,4,90,4).sum(axis=1).sum(axis=2)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

cmap = matplotlib.cm.get_cmap('viridis').copy()
cmap.set_bad(color = 'white', alpha = 0.2)

fig, [(ax1, ax2), (ax3, ax4)] = plt.subplots(2,2, figsize=(10, 5.5))

ax1 = plt.subplot(2,2,1,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im1 = ax1.imshow(cocco2/250/1000, extent=(-180,181,-90, 90), cmap=cmap, norm=LogNorm(vmin=0.01, vmax=10),  
                 origin='lower',  transform=ccrs.PlateCarree(), interpolation='nearest')
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax1.set_title('Coccolithophores', loc='left', fontsize=14)
ax1.set_ylabel('Latitude', fontsize=12)
ax1.set_yticks([-50,0,50])

ax2 = plt.subplot(2,2,2,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im2 = ax2.imshow(diatoms2/250/1000, extent=(-180,181,-90, 90), cmap=cmap, norm=LogNorm(vmin=0.01, vmax=10),  
                 origin='lower',  transform=ccrs.PlateCarree(), interpolation='nearest')
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax2.set_title('Diatoms', loc='left', fontsize=14)

ax3 = plt.subplot(2,2,3,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im3 = ax3.imshow(diazo2/250/1000, extent=(-180,181,-90, 90), cmap=cmap, norm=LogNorm(vmin=0.01, vmax=10),  
                 origin='lower',  transform=ccrs.PlateCarree(), interpolation='nearest')
ax3.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax3.set_yticks([-50,0,50])
ax3.set_xticks([-100,0,100])
ax3.set_ylabel('Latitude', fontsize=12)
ax3.set_xlabel('Longitude', fontsize=12)
ax3.set_title('Diazotrophs', loc='left', fontsize=14)

ax4 = plt.subplot(2,2,4,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im4 = ax4.imshow(pico2/250/1000, extent=(-180,181,-90, 90), cmap=cmap, norm=LogNorm(vmin=0.01, vmax=10),  
                 origin='lower',  transform=ccrs.PlateCarree(), interpolation='nearest')
ax4.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax4.set_xticks([-100,0,100])
ax4.set_xlabel('Longitude', fontsize=12)
ax4.set_title('Picophytoplankton', loc='left', fontsize=14)
fig.colorbar(im4,label='mg C $m^{-3}$', orientation="horizontal",  cax=fig.add_axes([0.4,-0.05,0.3,0.02]))

plt.tight_layout()
#plt.savefig('figures/FigureS8.pdf', bbox_inches='tight', transparent=True)
