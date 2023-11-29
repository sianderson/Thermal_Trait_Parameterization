#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Stephanie Anderson, Massachusetts Institute of Technology
Email: siander@mit.edu
    
This script examines estimates of export processes from each model simulation.

    INPUT: 
        Climate change simulations which have been archived here: https://doi.org/10.7910/DVN/6TLL8Z
        grid_igsm.nc: grid used in the model
        PFT Biomass with climate change (biomass_SQDUall.txt and biomass_DQDUall.txt)
    
    OUTPUT: 
        Figure S10. Mean biovolume
        Figure 5. Export processes
        
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

#%%  Extracting variables     
# In all files, SQSU = Eppley, SQDU = Kremer, DQDU = Anderson

# Climate change runs
SQDU_2080all_fname = 'Kremer_2080.nc'
SQDUall = nc.Dataset(SQDU_2080all_fname, 'r') 

DQDU_2080all_fname = 'Anderson_2080.nc'
DQDUall = nc.Dataset(DQDU_2080all_fname, 'r')

# Controls
SQDUc_fname = 'Kremer_control.nc'
SQDUc = nc.Dataset(SQDUc_fname, 'r')

DQDUc_fname = 'Anderson_control.nc'
DQDUc = nc.Dataset(DQDUc_fname, 'r')

# Particulate Organic Carbon (Export Production)
POC_SQDUall = SQDUall.variables['TRAC12']
POC_DQDUall = DQDUall.variables['TRAC12']

# Silica Export
Si_SQDUall = SQDUall.variables['TRAC16']
Si_DQDUall = DQDUall.variables['TRAC16']

# Particulate Inorganic Carbon
PIC_SQDUall = SQDUall.variables['TRAC17']
PIC_DQDUall = DQDUall.variables['TRAC17']


POC_SQDUc = SQDUc.variables['TRAC12']
POC_DQDUc = DQDUc.variables['TRAC12']

Si_SQDUc = SQDUc.variables['TRAC16']
Si_DQDUc = DQDUc.variables['TRAC16']

PIC_SQDUc = SQDUc.variables['TRAC17']
PIC_DQDUc = DQDUc.variables['TRAC17']

# Multiplying the particulate at 115m by the sinking rate
# Sinking rates are outlined in https://darwin3.readthedocs.io/
POC_sink_krem = POC_SQDUall[4] * 1.157407407407407E-004
POC_sink_and = POC_DQDUall[4] * 1.157407407407407E-004

POCc_sink_krem = POC_SQDUc[4] * 1.157407407407407E-004
POCc_sink_and = POC_DQDUc[4] * 1.157407407407407E-004 

Si_sink_krem = Si_SQDUall[4] * 1.157407407407407E-004
Si_sink_and = Si_DQDUall[4] * 1.157407407407407E-004

Sic_sink_krem = Si_SQDUc[4] * 1.157407407407407E-004
Sic_sink_and = Si_DQDUc[4] * 1.157407407407407E-004

Si_diff = np.subtract(Si_sink_and, Si_sink_krem)
Si_diffc = np.subtract(Sic_sink_and, Sic_sink_krem)

POC_diff = np.subtract(POC_sink_and, POC_sink_krem)
POC_diffc = np.subtract(POCc_sink_and, POCc_sink_krem)

PIC_diff = np.subtract(PIC_DQDUall[4], PIC_SQDUall[4])
PIC_diffc = np.subtract(PIC_DQDUc[4], PIC_SQDUc[4])


#%% Size Structure 

# Biovolumes for each phytoplankton phenotype in the model
BIOVOL  =  np.array([0.125892541179417     ,  0.410204102986607     ,   1.33659551654644     ,   4.35511873685569     ,
    14.1905752168909     ,   46.2381021399260     ,   150.660706618674     ,   490.907876152603     ,   1599.55802861467     ,
    14.1905752168909     ,   46.2381021399260     ,   150.660706618674     ,   490.907876152603     ,   1599.55802861467     ,
    14.1905752168909     ,   46.2381021399260     ,   150.660706618674     ,   490.907876152603     ,   1599.55802861467     ,
    5211.94711105080     ,   16982.4365246174     ,   55335.0109215737     ,   180301.774085957     ,   150.660706618674     ,
    490.907876152603     ,   1599.55802861467     ,   5211.94711105080     ,   16982.4365246174     ,   55335.0109215737     ,
    180301.774085957     ,   587489.352529777     ,   46.2381021399261     ,   150.660706618674     ,   490.907876152603     ,
    1599.55802861467     ,   5211.94711105081     ,   16982.4365246175     ,   55335.0109215737     ,   180301.774085957     ,
    587489.352529777     ,   1914255.92502109     ,   6237348.35482419     ,   20323570.1093622     ,   66221650.3701762     ,
    215774440.915267     ,   703072319.883834     ,   2290867652.76777     ,  3.863669770540691E-002,  0.125892541179417     ,
   0.410204102986607])
BIOVOL2 = BIOVOL[0:31]

ref2 = ['TRAC21','TRAC22','TRAC23','TRAC24','TRAC25','TRAC26','TRAC27','TRAC28','TRAC29',
        'TRAC30','TRAC31','TRAC32','TRAC33','TRAC34','TRAC35','TRAC36','TRAC37',
        'TRAC38','TRAC39','TRAC40','TRAC41','TRAC42','TRAC43',
       'TRAC44','TRAC45','TRAC46','TRAC47','TRAC48','TRAC49','TRAC50','TRAC51']

# biomass_SQDUall from q10_model_diversity_bray.py 
biomass_SQDUall = np.loadtxt('output/biomass_SQDUall.txt', dtype='float')
biomass_DQDUall = np.loadtxt('output/biomass_DQDUall.txt', dtype='float')

# for each PFT 
cell_biomass2_all = np.empty([len(ref2),90,144])
cell_biomass3_all = np.empty([len(ref2),90,144])

num = 0
np.seterr(invalid='warn')
for i in ref2: 

    cell2a = SQDUall.variables[i]
    cell3a = DQDUall.variables[i]
    cell_int2a = np.empty([6,90,144])
    cell_int3a = np.empty([6,90,144])
    
    for x in range(6): #240m
 
        cell_int2a[x] = depth_section[x] * cell2a[x]
        cell_int3a[x] = depth_section[x] * cell3a[x]
        
    cell_biomass2_all[num] = np.true_divide(np.nansum(cell_int2a, axis=0), biomass_SQDUall,
                                            out=np.zeros_like(np.nansum(cell_int2a, axis=0)), where=biomass_SQDUall!=0)
    cell_biomass3_all[num] = np.true_divide(np.nansum(cell_int3a, axis=0), biomass_DQDUall,
                                            out=np.zeros_like(np.nansum(cell_int3a, axis=0)), where=biomass_DQDUall!=0)
    num = num+1

cell_biomass2_all[cell_biomass2_all < 0.001]=np.nan
cell_biomass3_all[cell_biomass3_all < 0.001]=np.nan

#%% Biovolume (Figure S10)
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature

c2 = np.empty([31])
c3 = np.empty([31])
SQDU_vol = np.empty([90,144])
DQDU_vol = np.empty([90,144])

for j in range(90):
    for k in range(144):
        for i in range(31):
            c2[i] = cell_biomass2_all[i,j,k] * BIOVOL2[i]
            c3[i] = cell_biomass3_all[i,j,k] * BIOVOL2[i]
        SQDU_vol[j,k] = np.nansum(c2)/31
        DQDU_vol[j,k] = np.nansum(c3)/31

vol_diff = (DQDU_vol-SQDU_vol)/SQDU_vol*100

cmap2 = matplotlib.cm.get_cmap('viridis').copy()
cmap2.set_bad(color = '#440154FF', alpha = 1)

cmap3 = matplotlib.cm.get_cmap('RdBu_r').copy()
cmap3.set_bad(color = 'white', alpha = 1)

gridspec = dict(hspace=0.05)
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 7),
                                    subplot_kw={'projection': ccrs.Robinson(central_longitude=180)},
                                    gridspec_kw=gridspec)
im1 = ax1.imshow(SQDU_vol, extent=(0,360,-90,90), vmin=0, vmax=1000, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im1, label= 'Mean Biovolume \n($µm^3$)',ax=ax1, orientation='vertical', shrink=0.9, pad=0.02)

im2 = ax2.imshow(DQDU_vol, extent=(0,360,-90,90), vmin=0, vmax=1000, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im2, label= 'Mean Biovolume \n($µm^3$)',ax=ax2, orientation='vertical', shrink=0.9, pad=0.02)

im3 = ax3.imshow(vol_diff, extent=(0,360,-90,90), vmin=-100, vmax=100, origin='lower', cmap=cmap3, transform=ccrs.PlateCarree())
ax3.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im3, label= 'Difference \n(%)',ax=ax3, orientation='vertical', shrink=0.9, pad=0.02)

plt.gcf().text(0.47,0.73, 'Kremer', fontsize=12, weight='bold', rotation=90)
plt.gcf().text(0.47,0.45, 'Anderson', fontsize=12, weight='bold', rotation=90)
#plt.savefig('figures/FigureS10.pdf', bbox_inches='tight', transparent=True)


#%% Change in Export Production (Figure 5)

Si_a = (Si_sink_and - Sic_sink_and)
Si_k = (Si_sink_krem - Sic_sink_krem)

POC_a = (POC_sink_and - POCc_sink_and)
POC_k = (POC_sink_krem - POCc_sink_krem)

PIC_a = (PIC_DQDUall[4]-PIC_DQDUc[4])
PIC_k =  (PIC_SQDUall[4]-PIC_SQDUc[4])

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib

cmap2 = matplotlib.cm.get_cmap('RdBu_r').copy()
cmap2.set_bad(color = 'lightgrey', alpha = 1)

gridspec = dict(hspace=-0.15,height_ratios=[1,1,1],width_ratios=[1,1,1.2],wspace=0.05)
fig, ([ax1,ax4,ax7],[ax2,ax5,ax8],[ax3, ax6,ax9]) = plt.subplots(3, 3, figsize=(12, 7),
                                    subplot_kw={'projection': ccrs.Robinson(central_longitude=180)},
                                    gridspec_kw=gridspec)
# 1860
im1 = ax1.imshow(Si_diffc*0.028085*3.154e+7, extent=(0,360,-90,90), vmin=-10, vmax=10, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

im2 = ax2.imshow(POC_diffc*0.01201*3.154e+7, extent=(0,360,-90,90), vmin=-10, vmax=10, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))#fig.colorbar(im2, label= 'Export Production \n(gC/$m^{2}$/y)',ax=ax2, location="left",orientation='vertical', shrink=0.8, pad=0.02)

im3 = ax3.imshow(PIC_diffc*12.01, extent=(0,360,-90,90), vmin=-3, vmax=3, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax3.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m',facecolor='k'))#fig.colorbar(im3, label= 'PIC \n(mgC/$m^{3}$)',ax=ax3, location="left",orientation='vertical', shrink=0.8, pad=0.02)


# 2100
## Kremer
im4 = ax4.imshow((Si_k)*0.028085*3.154e+7, extent=(0,360,-90,90), vmin=-10, vmax=10, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax4.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

im5 = ax5.imshow((POC_k)*0.01201*3.154e+7, extent=(0,360,-90,90), vmin=-10, vmax=10, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax5.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

im6 = ax6.imshow(PIC_k*12.01, extent=(0,360,-90,90), vmin=-3, vmax=3, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax6.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m',facecolor='k'))

## Anderson
im7 = ax7.imshow(Si_a*0.028085*3.154e+7, extent=(0,360,-90,90), vmin=-10, vmax=10, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax7.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im7, label= 'Silica Export \n(g Si $m^{-2}$ $y^{-1}$)',ax=ax7, orientation='vertical', shrink=0.7, pad=0.02)

im8 = ax8.imshow(POC_a*0.01201*3.154e+7, extent=(0,360,-90,90), vmin=-10, vmax=10, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax8.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im8, label= 'Export Production \n(g C $m^{-2}$ $y^{-1}$)',ax=ax8, orientation='vertical', shrink=0.7, pad=0.02)

im9 = ax9.imshow(PIC_a*12.01, extent=(0,360,-90,90), vmin=-3, vmax=3, origin='lower', cmap=cmap2, transform=ccrs.PlateCarree())
ax9.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m',facecolor='k'))
fig.colorbar(im9, label= 'PIC \n(mg C $m^{-3}$)',ax=ax9, orientation='vertical', shrink=0.7, pad=0.02)

line=plt.Line2D((0.365,0.365),(0.15,0.9), color='grey', linestyle='dashed')
fig.add_artist(line)

plt.gcf().text(0.12,0.83, '(a)', fontsize=12, weight='bold')
plt.gcf().text(0.37,0.83, '(b)', fontsize=12, weight='bold')
plt.gcf().text(0.62,0.83, '(c)', fontsize=12, weight='bold')

plt.gcf().text(0.12,0.59, '(d)', fontsize=12, weight='bold')
plt.gcf().text(0.37,0.59, '(e)', fontsize=12, weight='bold')
plt.gcf().text(0.62,0.59, '(f)', fontsize=12, weight='bold')

plt.gcf().text(0.12,0.35, '(g)', fontsize=12, weight='bold')
plt.gcf().text(0.37,0.35, '(h)', fontsize=12, weight='bold')
plt.gcf().text(0.62,0.35, '(i)', fontsize=12, weight='bold')

plt.gcf().text(0.22,0.9, '1860', fontsize=14, weight='bold')
plt.gcf().text(0.18,0.87, '(Anderson-Kremer)', fontsize=12)
plt.gcf().text(0.45,0.9, '$\Delta$ Kremer', fontsize=14, weight='bold')
plt.gcf().text(0.45,0.87, '(2100-1860)', fontsize=12)
plt.gcf().text(0.68,0.9, '$\Delta$ Anderson', fontsize=14, weight='bold')
plt.gcf().text(0.69,0.87, '(2100-1860)', fontsize=12)
#plt.savefig('figures/Figure5.pdf', bbox_inches='tight', transparent=True)
