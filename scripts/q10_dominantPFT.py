#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Stephanie Anderson, Massachusetts Institute of Technology
Email: siander@mit.edu
    
This script examines the dominant PFT globally, and the PFT latitudinal extent differences between each model.

    INPUT: 
        Control simulations which have been archived here: https://doi.org/10.7910/DVN/6TLL8Z
        grid_igsm.nc: grid used in the model
    
    OUTPUT: 
        Figure 3C. Global maps of the dominant PFT for each model.
        
%=========================================================================
"""

#%%
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# depth calculations
ds = nc.Dataset('data/grid_igsm.nc', 'r') 
depth = ds.variables['Z'][:]
depth = depth*-1

# Calculate depth of each section
depth2 = np.empty([22])
depth_section = np.empty([22])
for x in range(22):
    if x==0:
        depth2[x] = 2*depth[x]
        depth_section[x] = depth2[x]
    if x>0:
        depth2[x] = 2*depth[x]-depth2[x-1]
        depth_section[x] = depth2[x] - depth2[x-1]

#%% Depth integrated primary production (PP)
# Controls
# In all files, SQSU = Eppley, SQDU = Kremer, DQDU = Anderson
SQSUc_fname = 'Eppley_control.nc'
SQSUc = nc.Dataset(SQSUc_fname, 'r')

SQDUc_fname = 'Kremer_control.nc'
SQDUc = nc.Dataset(SQDUc_fname, 'r')

DQDUc_fname = 'Anderson_control.nc'
DQDUc = nc.Dataset(DQDUc_fname, 'r')

# creating a 'land' mask to account for NAs in model output
land_fname = 'Anderson_control.nc'
land = nc.Dataset(land_fname, 'r')

# get primary production
land = land.variables['PP']
PP_SQSUc = SQSUc.variables['PP']
PP_SQDUc = SQDUc.variables['PP']
PP_DQDUc = DQDUc.variables['PP']

PP_int_SQSUc = np.empty([22,90,144])
PP_int_sum_SQSUc = np.empty([22,90,144])

PP_int_SQDUc = np.empty([22,90,144])
PP_int_sum_SQDUc = np.empty([22,90,144])
PP_int_DQDUc = np.empty([22,90,144])
PP_int_sum_DQDUc = np.empty([22,90,144])

for x in range(22):
    PP_int_SQSUc[x] = depth_section[x] * PP_SQSUc[x]
    PP_int_SQDUc[x] = depth_section[x] * PP_SQDUc[x]
    PP_int_DQDUc[x] = depth_section[x] * PP_DQDUc[x]
    if x<1:
        PP_int_sum_SQSUc[x]=PP_int_SQSUc[x]
        PP_int_sum_SQDUc[x]=PP_int_SQDUc[x]
        PP_int_sum_DQDUc[x]=PP_int_DQDUc[x]
    if x>=1:
        PP_int_sum_SQSUc[x]=PP_int_sum_SQSUc[x-1]+PP_int_SQSUc[x]
        PP_int_sum_SQDUc[x]=PP_int_sum_SQDUc[x-1]+PP_int_SQDUc[x]
        PP_int_sum_DQDUc[x]=PP_int_sum_DQDUc[x-1]+PP_int_DQDUc[x]
     
#%% Average PP over top 240 m
ref = {'coccolithophores':['TRAC25','TRAC26','TRAC27','TRAC28','TRAC29'],
       'cyano':['TRAC21','TRAC22'],
       'diatoms':['TRAC43','TRAC42','TRAC41','TRAC40','TRAC39','TRAC38','TRAC37','TRAC36','TRAC35'],
       'diazotrophs':['TRAC30','TRAC31','TRAC32','TRAC33','TRAC34'],
       'dinos':['TRAC44','TRAC45','TRAC46','TRAC47','TRAC48','TRAC49','TRAC50','TRAC51'],
       'greens':['TRAC23','TRAC24'],
       }

pft_sum_SQSU = np.empty([len(ref),90,144])
pft_sum_SQDU = np.empty([len(ref),90,144])
pft_sum_DQDU = np.empty([len(ref),90,144])
num1 = 0

for key in ref: # for each PFT
    all_cell_sum = np.empty([len(ref[key]),90,144])
    all_cell_sum2 = np.empty([len(ref[key]),90,144])
    all_cell_sum3 = np.empty([len(ref[key]),90,144])
    num2 = 0
    for i in ref[key]: # for each phenotype (tracer)
        cell = SQSUc.variables[i]
        cell2 = SQDUc.variables[i]
        cell3 = DQDUc.variables[i]
        cell_int = np.empty([8,90,144])   
        cell_sum = np.empty([9,90,144])
        cell_int2 = np.empty([8,90,144])   
        cell_sum2 = np.empty([9,90,144])
        cell_int3 = np.empty([8,90,144])   
        cell_sum3 = np.empty([9,90,144])
        for x in range(6): # 6 to 240 m
            cell_int[x] = depth_section[x] * cell[x]
            cell_int2[x] = depth_section[x] * cell2[x]
            cell_int3[x] = depth_section[x] * cell3[x]
            if x<1:
                cell_sum[x]=cell_int[x]
                cell_sum2[x]=cell_int2[x]
                cell_sum3[x]=cell_int3[x]
            if x>=1:
                cell_sum[x]=cell_sum[x-1]+cell_int[x]
                cell_sum2[x]=cell_sum2[x-1]+cell_int2[x]
                cell_sum3[x]=cell_sum3[x-1]+cell_int3[x]
        all_cell_sum[num2] = cell_sum[3]
        all_cell_sum2[num2] = cell_sum2[3]
        all_cell_sum3[num2] = cell_sum3[3]
        num2 = num2+1
    pft_sum_SQSU[num1] = np.sum(all_cell_sum, axis=0)
    pft_sum_SQDU[num1] = np.sum(all_cell_sum2, axis=0)
    pft_sum_DQDU[num1] = np.sum(all_cell_sum3, axis=0)
    num1 = num1+1

biomass_DQDU = np.nansum(pft_sum_DQDU, axis=0)
biomass_SQDU = np.nansum(pft_sum_SQDU, axis=0)
biomass_SQSU = np.nansum(pft_sum_SQSU, axis=0)


#%% Dominant PFT (Figure 3C)

from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

pft_max_DQDU = pft_sum_DQDU.argmax(axis=0)
pft_max_SQSU = pft_sum_SQSU.argmax(axis=0)
pft_max_SQDU = pft_sum_SQDU.argmax(axis=0)

SQSU_land = np.ma.masked_where(biomass_SQSU < 0.001, pft_max_SQSU)
SQDU_land = np.ma.masked_where(biomass_SQDU < 0.001, pft_max_SQDU)
DQDU_land = np.ma.masked_where(biomass_DQDU < 0.001, pft_max_DQDU)

pfts = (["Coccolithophores", "Cyanobacteria", "Diatoms","Diazotrophs","Dinoflagellates", "Green Algae"])
colors2=["orange","#ec3a25","#026cb1", "purple", "brown", "#3ea127",  "green"]
nodes = [0.0, 0.17, 0.35, 0.4,0.6,0.85,1.0]
cmap1 = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes,colors2)))
cmap1.set_bad('black', 0.5)

gridspec = dict(hspace=0.2,bottom=0, height_ratios=[1, 1,1])
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, 
                                    figsize=(12, 7), 
                                    gridspec_kw=gridspec,sharex='col',
                                    subplot_kw={'projection': ccrs.Robinson(central_longitude=180)})

im1 = ax1.imshow(SQSU_land, extent=(0,360,-90,90), vmin=0, vmax=6, cmap=cmap1,
                 origin='lower', transform=ccrs.PlateCarree())
ax1.set_title('Eppley', loc='left', fontweight='bold')
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

im2 = ax2.imshow(SQDU_land, extent=(0,360,-90,90), vmin=0, vmax=6, cmap=cmap1,
                 origin='lower', transform=ccrs.PlateCarree())
ax2.set_ylabel('Latitude', fontsize=12)
ax2.set_title('Kremer', loc='left', fontweight='bold')
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

im3 = ax3.imshow(DQDU_land, extent=(0,360,-90,90), vmin=0, vmax=6, cmap=cmap1,
                 origin='lower', transform=ccrs.PlateCarree())
ax3.set_xlabel('Longitude', fontsize=12)
ax3.set_title('Anderson', loc='left', fontweight='bold')
ax3.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

values = range(6)
colors = [im1.cmap(im1.norm(value)) for value in values]

# create a patch for every color 
patches = [mpatches.Patch(color=colors[i], label=pfts[i]) for i in range(len(values))]
plt.legend(handles=patches, bbox_to_anchor=(1.3, 1.4), loc="lower center", borderaxespad=0., ncol=1, frameon=False)
plt.tight_layout()

#plt.savefig('figures/Figure3C.pdf', bbox_inches='tight', transparent=True)

#%% Extent changes

import math
# Difference in area covered
xg = ds.variables['XG'][:] #longitudinal coordinate at point
yg = ds.variables['YG'][:]  #latitudinal coordinate at point

# input is in degrees
# output is in km^2
def surface(lat1, lon1, lat2, lon2):
    rad = math.pi/180
    A1 = 2 * math.pi * (1-math.sin(lat1*rad)) * (6371 ** 2)
    A2 = 2 * math.pi * (1-math.sin(lat2*rad)) * (6371 ** 2)
    surf = abs(A1-A2)*abs(lon1 - lon2)/360
    return(surf)

surfacearea = np.empty([90,144])
for i in range(90):
    lat1 = yg[i+1,0]
    lat2 = yg[i,0]
    lon1 = xg[0,1]
    lon2 = xg[0,0]
    sa = surface(lat1, lon1, lat2, lon2)
    surfacearea[i,:] = sa

# Coccolithophores
sqdu_bi = np.where(pft_sum_SQDU[0] > 0.001, 1, 0) 
dqdu_bi = np.where(pft_sum_DQDU[0] > 0.001, 1, 0) 
per = np.where((dqdu_bi - sqdu_bi) <= 0, 0, 1) 
(per*surfacearea).sum(axis=None)/(sqdu_bi*surfacearea).sum(axis=None)

# Cyanobacteria
sqdu_bi = np.where(pft_sum_SQDU[1] > 0.001, 1, 0) 
dqdu_bi = np.where(pft_sum_DQDU[1] > 0.001, 1, 0) 
per = np.where((sqdu_bi - dqdu_bi) <= 0, 0, 1) 
(per*surfacearea).sum(axis=None)/(sqdu_bi*surfacearea).sum(axis=None)
