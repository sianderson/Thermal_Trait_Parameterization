#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Stephanie Anderson, Massachusetts Institute of Technology
Email: siander@mit.edu
    
This script examines community change using Bray-Curtis dissimilarity.

    INPUT: 
        Control simulations which has been archived here ############
        grid_igsm.nc: grid used in the model
    
    OUTPUT: 
       Figure 6. Bray-Curtis dissimilarity
       Figure 7. Mean change in depth-averaged PFT biomass
       Figure S6. Mean change in each PFT biomass
       PFT biomass with climate change (biomass_SQDUall.txt and biomass_DQDUall.txt)
        
%=========================================================================
"""
#%%
import os 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

os.chdir("/Users/Stephanie/Desktop/MIT/Q10_Variability/Code & Datasets/")
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


SQDUc_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_SQDU_control_avg.nc'
SQDUc = nc.Dataset(SQDUc_fname, 'r')
DQDUc_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_control_avg.nc'
DQDUc = nc.Dataset(DQDUc_fname, 'r')

SQDU_2080all_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_SQDU_2080avg.nc'
SQDUall = nc.Dataset(SQDU_2080all_fname, 'r') 
DQDU_2080all_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_2080avg.nc'
DQDUall = nc.Dataset(DQDU_2080all_fname, 'r')

PP_SQDUc = SQDUc.variables['PP']
PP_DQDUc = DQDUc.variables['PP']
PP_SQDUall = SQDUall.variables['PP']
PP_DQDUall = DQDUall.variables['PP']

#%% Depth integrated PP

# Controls
PP_int_SQDUc = np.empty([22,90,144])
PP_int_sum_SQDUc = np.empty([22,90,144])
PP_int_DQDUc = np.empty([22,90,144])
PP_int_sum_DQDUc = np.empty([22,90,144])

for x in range(22):
    PP_int_SQDUc[x] = depth_section[x] * PP_SQDUc[x]
    PP_int_DQDUc[x] = depth_section[x] * PP_DQDUc[x]
    if x<1:
        PP_int_sum_SQDUc[x]=PP_int_SQDUc[x]
        PP_int_sum_DQDUc[x]=PP_int_DQDUc[x]
    if x>=1:
        PP_int_sum_SQDUc[x]=PP_int_sum_SQDUc[x-1]+PP_int_SQDUc[x]
        PP_int_sum_DQDUc[x]=PP_int_sum_DQDUc[x-1]+PP_int_DQDUc[x]

# Climate Change Scenarios
PP_int_SQDUall = np.empty([22,90,144])
PP_int_sum_SQDUall = np.empty([22,90,144])
PP_int_DQDUall = np.empty([22,90,144])
PP_int_sum_DQDUall = np.empty([22,90,144])

for x in range(22):
    PP_int_SQDUall[x] = depth_section[x] * PP_SQDUall[x]
    PP_int_DQDUall[x] = depth_section[x] * PP_DQDUall[x]
    if x<1:
        PP_int_sum_SQDUall[x]=PP_int_SQDUall[x]
        PP_int_sum_DQDUall[x]=PP_int_DQDUall[x]
    if x>=1:
        PP_int_sum_SQDUall[x]=PP_int_sum_SQDUall[x-1]+PP_int_SQDUall[x]
        PP_int_sum_DQDUall[x]=PP_int_sum_DQDUall[x-1]+PP_int_DQDUall[x]
  
#%% Average biomass over top 260 m

ref = {'coccolithophores':['TRAC25','TRAC26','TRAC27','TRAC28','TRAC29'],
       'cyano':['TRAC21','TRAC22'],
       'diatoms':['TRAC43','TRAC42','TRAC41','TRAC40','TRAC39','TRAC38','TRAC37','TRAC36','TRAC35'],
       'diazotrophs':['TRAC30','TRAC31','TRAC32','TRAC33','TRAC34'],
       'dinos':['TRAC44','TRAC45','TRAC46','TRAC47','TRAC48','TRAC49','TRAC50','TRAC51'],
       'greens':['TRAC23','TRAC24'],
       }

pft_sum_SQDUa = np.empty([len(ref),90,144])
pft_sum_DQDUa = np.empty([len(ref),90,144])
num1 = 0

for key in ref: # for each PFT
    all_cell_sum2a = np.empty([len(ref[key]),90,144])
    all_cell_sum3a = np.empty([len(ref[key]),90,144])
    num2 = 0
    for i in ref[key]: # for each cell (tracer)
        cell2a = SQDUall.variables[i]
        cell3a = DQDUall.variables[i]
        cell_int2a = np.empty([8,90,144])   
        cell_sum2a = np.empty([9,90,144])
        cell_int3a = np.empty([8,90,144])   
        cell_sum3a = np.empty([9,90,144])
        for x in range(6): # 6 to 240m
            cell_int2a[x] = depth_section[x] * cell2a[x]
            cell_int3a[x] = depth_section[x] * cell3a[x]
            if x<1:
                cell_sum2a[x]=cell_int2a[x]
                cell_sum3a[x]=cell_int3a[x]
            if x>=1:
                cell_sum2a[x]=cell_sum2a[x-1]+cell_int2a[x]
                cell_sum3a[x]=cell_sum3a[x-1]+cell_int3a[x]
        all_cell_sum2a[num2] = cell_sum2a[3]
        all_cell_sum3a[num2] = cell_sum3a[3]
        num2 = num2+1
    pft_sum_SQDUa[num1] = np.sum(all_cell_sum2a, axis=0)
    pft_sum_DQDUa[num1] = np.sum(all_cell_sum3a, axis=0)
    num1 = num1+1

biomass_DQDUall = np.nansum(pft_sum_DQDUa, axis=0)
biomass_SQDUall = np.nansum(pft_sum_SQDUa, axis=0)

pft_biomass_SQDUa = np.empty([len(ref),90,144])
pft_biomass_DQDUa = np.empty([len(ref),90,144])
for i in range(6):
    pft_biomass_SQDUa[i] = np.divide(pft_sum_SQDUa[i],biomass_SQDUall)
    pft_biomass_DQDUa[i] = np.divide(pft_sum_DQDUa[i],biomass_DQDUall)
    
#np.savetxt('output/biomass_SQDUall.txt', biomass_SQDUall, fmt='%.4e')
#np.savetxt('output/biomass_DQDUall.txt', biomass_DQDUall, fmt='%.4e')


#%% Average over top 260 m
ref = {'coccolithophores':['TRAC25','TRAC26','TRAC27','TRAC28','TRAC29'],
       'cyano':['TRAC21','TRAC22'],
       'diatoms':['TRAC43','TRAC42','TRAC41','TRAC40','TRAC39','TRAC38','TRAC37','TRAC36','TRAC35'],
       'diazotrophs':['TRAC30','TRAC31','TRAC32','TRAC33','TRAC34'],
       'dinos':['TRAC44','TRAC45','TRAC46','TRAC47','TRAC48','TRAC49','TRAC50','TRAC51'],
       'greens':['TRAC23','TRAC24'],
       }

pft_sum_SQDU = np.empty([len(ref),90,144])
pft_sum_DQDU = np.empty([len(ref),90,144])
num1 = 0

for key in ref: # for each PFT
    all_cell_sum2 = np.empty([len(ref[key]),90,144])
    all_cell_sum3 = np.empty([len(ref[key]),90,144])
    num2 = 0
    for i in ref[key]: # for each phenotype (tracer)
        cell2 = SQDUc.variables[i]
        cell3 = DQDUc.variables[i]
        cell_int2 = np.empty([8,90,144])   
        cell_sum2 = np.empty([9,90,144])
        cell_int3 = np.empty([8,90,144])   
        cell_sum3 = np.empty([9,90,144])
        for x in range(6): # 6 to 240m
            cell_int2[x] = depth_section[x] * cell2[x]
            cell_int3[x] = depth_section[x] * cell3[x]
            if x<1:
                cell_sum2[x]=cell_int2[x]
                cell_sum3[x]=cell_int3[x]
            if x>=1:
                cell_sum2[x]=cell_sum2[x-1]+cell_int2[x]
                cell_sum3[x]=cell_sum3[x-1]+cell_int3[x]
        all_cell_sum2[num2] = cell_sum2[3]
        all_cell_sum3[num2] = cell_sum3[3]
        num2 = num2+1
    pft_sum_SQDU[num1] = np.sum(all_cell_sum2, axis=0)
    pft_sum_DQDU[num1] = np.sum(all_cell_sum3, axis=0)
    num1 = num1+1

biomass_DQDU = np.nansum(pft_sum_DQDU, axis=0)
biomass_SQDU = np.nansum(pft_sum_SQDU, axis=0)

biomass_SQDU[biomass_SQDU < 0.001]=np.nan
biomass_DQDU[biomass_DQDU < 0.001]=np.nan

pft_biomass_SQDU = np.empty([len(ref),90,144])
pft_biomass_DQDU = np.empty([len(ref),90,144])
for i in range(6):
    pft_biomass_SQDU[i] = np.divide(pft_sum_SQDU[i],biomass_SQDU)
    pft_biomass_DQDU[i] = np.divide(pft_sum_DQDU[i],biomass_DQDU)
    
#%% Change in community structure (PFT)
# bray curtis calculations
        
c2 = np.empty([6])
d2 = np.empty([6])
c3 = np.empty([6])
d3 = np.empty([6])

c4 = np.empty([6])
d4 = np.empty([6])
c5 = np.empty([6])
d5 = np.empty([6])
SQDU_bray = np.empty([90,144])
DQDU_bray = np.empty([90,144])
bray_control = np.empty([90,144])
bray_all = np.empty([90,144])

for j in range(90):
    for k in range(144):
        for i in range(6):
            c2[i]=abs(pft_biomass_SQDUa[i,j,k] - pft_biomass_SQDU[i,j,k])
            d2[i]=pft_biomass_SQDUa[i,j,k] + pft_biomass_SQDU[i,j,k]
            c3[i]=abs(pft_biomass_DQDUa[i,j,k] - pft_biomass_DQDU[i,j,k])
            d3[i]=pft_biomass_DQDUa[i,j,k] + pft_biomass_DQDU[i,j,k]
            
            c4[i]=abs(pft_biomass_DQDU[i,j,k] - pft_biomass_SQDU[i,j,k])
            d4[i]=pft_biomass_DQDU[i,j,k] + pft_biomass_SQDU[i,j,k]
            c5[i]=abs(pft_biomass_DQDUa[i,j,k] - pft_biomass_SQDUa[i,j,k])
            d5[i]=pft_biomass_DQDUa[i,j,k] + pft_biomass_SQDUa[i,j,k]
        SQDU_bray[j,k]=np.divide(np.nansum(c2),np.nansum(d2),
                                    out=np.zeros_like(np.nansum(c2)), where=np.nansum(d2)!=0)
        DQDU_bray[j,k]=np.divide(np.nansum(c3),np.nansum(d3),
                                     out=np.zeros_like(np.nansum(c3)), where=np.nansum(d3)!=0)
        
        bray_control[j,k]=np.divide(np.nansum(c4),np.nansum(d4),
                                    out=np.zeros_like(np.nansum(c4)), where=np.nansum(d4)!=0)
        bray_all[j,k]=np.divide(np.nansum(c5),np.nansum(d5),
                                     out=np.zeros_like(np.nansum(c5)), where=np.nansum(d5)!=0)


#%% Bray curtis (Figure 6)
import matplotlib
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

land_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_control_avg.nc'
land = nc.Dataset(land_fname, 'r')
land = land.variables['PP']

cmap3 = matplotlib.cm.get_cmap('RdBu').copy()
cmap3.set_bad(color = 'white', alpha = 1)
vmax=0.4

SQDU_braym =np.ma.masked_where((land[0] == 0), SQDU_bray)
DQDU_braym =np.ma.masked_where((land[0] == 0), DQDU_bray)
bray_controlm =np.ma.masked_where((land[0] == 0), bray_control)
bray_allm =np.ma.masked_where((land[0] == 0), bray_all)

# figure
bray_med_SQDU = np.nanmean(SQDU_braym.filled(np.nan), axis=1)
bray_std_SQDU = np.nanstd(SQDU_braym.filled(np.nan), axis=1) 
bray_med_DQDU = np.nanmean(DQDU_braym.filled(np.nan), axis=1)
bray_std_DQDU = np.nanstd(DQDU_braym.filled(np.nan), axis=1) 
bray_med_con = np.nanmean(bray_controlm.filled(np.nan), axis=1)
bray_std_con = np.nanstd(bray_controlm.filled(np.nan), axis=1) 
bray_med_all = np.nanmean(bray_allm.filled(np.nan), axis=1)
bray_std_all = np.nanstd(bray_allm.filled(np.nan), axis=1) 

cmap2 = matplotlib.cm.get_cmap('cividis').copy()
cmap2.set_bad(color = 'k', alpha = 1)
vmax=0.4

gridspec = dict(bottom=0.15,height_ratios=[1,1],wspace=0.1, width_ratios=[1,1,0.5])
fig, [(ax5,ax7, ax9),(ax6, ax8, ax9)] = plt.subplots(2,3, figsize=(9,4), sharey='row', gridspec_kw=gridspec)

ax5 = plt.subplot2grid((2, 3), (0, 0), projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im5= ax5.imshow(SQDU_bray, extent=(0,360,-90,90), cmap=cmap2, vmin=0, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax5.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax5.set_yticks([-50,0,50])
ax5.set_title('Kremer',loc='left')

ax6 = plt.subplot2grid((2, 3), (1, 0), projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im6= ax6.imshow(DQDU_bray, extent=(0,360,-90,90), cmap=cmap2, vmin=0, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax6.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax6.set_title('Anderson',loc='left')
ax6.set_xticks([-100,0,100])
ax6.set_xticklabels([80,180,280])
ax6.set_yticks([-50,0,50])

ax7 = plt.subplot2grid((2, 3), (0, 1), projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im7= ax7.imshow(bray_control, extent=(0,360,-90,90), cmap=cmap2, vmin=0, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax7.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax7.set_title('Control',loc='left')


ax8 = plt.subplot2grid((2, 3), (1, 1), projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im8= ax8.imshow(bray_all, extent=(0,360,-90,90), cmap=cmap2, vmin=0, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax8.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im8,label='Bray-Curtis Dissimilarity', shrink=0.65, 
             orientation='horizontal', norm=colors.CenteredNorm(),cax=fig.add_axes([0.29,-0.08,0.3,0.03]))
ax8.set_title('Climate Change',loc='left')
ax8.set_xticks([-100,0,100])
ax8.set_xticklabels([80,180,280])


c2 = '#2A788EFF'
c3 = '#7AD151FF'

ax1 = plt.subplot2grid((2, 3), (0, 2), rowspan=2)
ax1.plot(bray_med_SQDU,range(-90,90,2), color=c2)
ax1.fill_betweenx(range(-90,90,2),bray_med_SQDU-bray_std_SQDU, bray_med_SQDU+bray_std_SQDU, facecolor=c2, alpha=0.1)
ax1.fill_betweenx(range(-90,90,2),bray_med_DQDU-bray_std_DQDU, bray_med_DQDU+bray_std_DQDU, facecolor=c3, alpha=0.1)
ax1.fill_betweenx(range(-90,90,2),bray_med_con-bray_std_con, bray_med_con+bray_std_con, facecolor='k', alpha=0.1)
ax1.fill_betweenx(range(-90,90,2),bray_med_all-bray_std_all, bray_med_all+bray_std_all, facecolor='red', alpha=0.1)
ax1.axhline(y=0, color='grey', linestyle='--')
line1, = ax1.plot(bray_med_SQDU,range(-90,90,2), color=c2, label='Kremer')
line2, = ax1.plot(bray_med_DQDU,range(-90,90,2), color=c3, label='Anderson')
line3, = ax1.plot(bray_med_con,range(-90,90,2), color='k', linestyle='dashed', label='Control')
line4, = ax1.plot(bray_med_all,range(-90,90,2), color='red',linestyle='dashed', label='2100all')
ax1.set_xlim(-0.1,0.6)
ax1.set_yticks([0])
leg = ax1.legend(handles = [line1,line2,line3, line4],bbox_to_anchor=(0.55, -0.5), loc="lower center", frameon=False)
for line in leg.get_lines():
    line.set_linewidth(3.0)
ax1.set_xlabel('Bray Curtis', fontsize=11)

fig.text(0.07,0.5, 'Latitude', va='center', rotation='vertical', size=12)
fig.text(0.4,0.04, 'Longitude', va='center', size=12)
plt.gcf().text(0.125,0.95, '(a) 1860 vs 2100', fontsize=12, weight='bold')
plt.gcf().text(0.44,0.95, '(b) Anderson vs Kremer', fontsize=12, weight='bold')
plt.gcf().text(0.755,0.95, '(c)', fontsize=12, weight='bold')

#plt.savefig('figures/Figure6.pdf', bbox_inches='tight', transparent=True)

#%%# PFT Biomass change 2100
# Figure S6

import cartopy.crs as ccrs
import cartopy.feature as cfeature

d=240 # depth
vmax=1
vmin=-1

cmap2 = matplotlib.cm.get_cmap('RdBu_r').copy()
cmap2.set_bad(color = 'lightgrey', alpha = 1)

gridspec = dict(hspace=0.1,bottom=0.05, width_ratios=[1,1,1], wspace=-0.7)
fig, ([ax1, ax7, ax13], [ax2, ax8,ax14], [ax3, ax9, ax15],
      [ax4,ax10,ax16], [ax5,ax11, ax17], [ax6,ax12,ax18]) = plt.subplots(6,3, figsize=(11, 9), sharey='col', sharex='col',gridspec_kw=gridspec)
                                                                                        

ax1 = plt.subplot(641,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im1 = ax1.imshow(((pft_sum_SQDUa[0]-pft_sum_SQDU[0])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                 origin='lower', transform=ccrs.PlateCarree())
ax1.set_yticks([-50,0,50])
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax1.set_facecolor('black')
ax1.set_title('Coccolithophores', loc='left')


ax2 = plt.subplot(645,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im2= ax2.imshow(((pft_sum_SQDUa[1]-pft_sum_SQDU[1])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax2.set_yticks([-50,0,50])
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax2.set_facecolor('black')
ax2.set_title('Cyanobacteria', loc='left')

ax3 = plt.subplot(649,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im3= ax3.imshow(((pft_sum_SQDUa[2]-pft_sum_SQDU[2])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax3.set_yticks([-50,0,50])
ax3.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax3.set_facecolor('black')
ax3.set_title('Diatoms', loc='left')

ax4 = plt.subplot(6,4,13,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im4= ax4.imshow(((pft_sum_SQDUa[3]-pft_sum_SQDU[3])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax4.set_yticks([-50,0,50])
ax4.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax4.set_ylabel('Latitude', fontsize=12)
ax4.yaxis.set_label_coords(-0.15,1.07)
ax4.set_facecolor('black')
ax4.set_yticks([-50,0,50])
ax4.set_yticklabels([50,0,-50])
ax4.set_title('Diazotrophs', loc='left')

ax5 = plt.subplot(6,4,17,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im5= ax5.imshow(((pft_sum_SQDUa[4]-pft_sum_SQDU[4])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax5.set_yticks([-50,0,50])
ax5.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax5.set_facecolor('black')
ax5.set_title('Dinoflagellates', loc='left')

ax6 = plt.subplot(6,4,21,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im6= ax6.imshow(((pft_sum_SQDUa[5]-pft_sum_SQDU[5])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax6.set_yticks([-50,0,50])
ax6.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax6.set_facecolor('black')
ax6.set_title('Green Algae', loc='left')
ax6.set_xticks([-100,0,100])
ax6.set_xticklabels([80,180,280])

ax7 = plt.subplot(6,4,2,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im7 = ax7.imshow(((pft_sum_DQDUa[0]-pft_sum_DQDU[0])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax7.set_yticks([])
ax7.set_facecolor('black')
ax7.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

ax8 = plt.subplot(6,4,6,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im8=  ax8.imshow(((pft_sum_DQDUa[1]-pft_sum_DQDU[1])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax8.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax8.set_yticks([])
ax8.set_facecolor('black')

ax9 = plt.subplot(6,4,10,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im9=  ax9.imshow(((pft_sum_DQDUa[2]-pft_sum_DQDU[2])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax9.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax9.set_yticks([])
ax9.set_facecolor('black')

ax10 = plt.subplot(6,4,14,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im10=  ax10.imshow(((pft_sum_DQDUa[3]-pft_sum_DQDU[3])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax10.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax10.set_facecolor('black')
ax10.set_yticks([])

ax11 = plt.subplot(6,4,18,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im11=  ax11.imshow(((pft_sum_DQDUa[4]-pft_sum_DQDU[4])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax11.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax11.set_yticks([])
ax11.set_facecolor('black')

ax12 = plt.subplot(6,4,22,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im12=  ax12.imshow(((pft_sum_DQDUa[5]-pft_sum_DQDU[5])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                origin='lower', transform=ccrs.PlateCarree())
ax12.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax12.set_yticks([])
ax12.set_xlabel('Longitude', fontsize=12)
ax12.set_xticks([-100,0,100])
ax12.set_xticklabels([80,180,280])

ax13 = plt.subplot(6,4,3,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im13 = ax13.imshow(((pft_sum_DQDUa[0]-pft_sum_SQDUa[0])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                   origin='lower', transform=ccrs.PlateCarree())#norm=matplotlib.colors.LogNorm())
ax13.set_yticks([])
ax13.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax13.set_facecolor('black')

ax14 = plt.subplot(6,4,7,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im14= ax14.imshow(((pft_sum_DQDUa[1]-pft_sum_SQDUa[1])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                  origin='lower', transform=ccrs.PlateCarree())
ax14.set_yticks([])
ax14.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax14.set_yticks([])
ax14.set_facecolor('black')

ax15 = plt.subplot(6,4,11,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im15= ax15.imshow(((pft_sum_DQDUa[2]-pft_sum_SQDUa[2])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                  origin='lower', transform=ccrs.PlateCarree())
ax15.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax15.set_yticks([])
ax15.set_facecolor('black')

ax16 = plt.subplot(6,4,15,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im16= ax16.imshow(((pft_sum_DQDUa[3]-pft_sum_SQDUa[3])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                  origin='lower', transform=ccrs.PlateCarree())
ax16.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax16.set_facecolor('black')
ax16.set_yticks([])

ax17 = plt.subplot(6,4,19,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im17= ax17.imshow(((pft_sum_DQDUa[4]-pft_sum_SQDUa[4])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                  origin='lower', transform=ccrs.PlateCarree())
ax17.set_yticks([])
ax17.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax17.set_facecolor('black')

ax18 = plt.subplot(6,4,23,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im18= ax18.imshow(((pft_sum_DQDUa[5]-pft_sum_SQDUa[5])*12.01/d), extent=(0,360,-90,90), vmin=vmin, vmax=vmax, cmap=cmap2,
                  origin='lower', transform=ccrs.PlateCarree())
ax18.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im1,label='mg C $m^{-3}$', orientation="horizontal",  cax=fig.add_axes([0.27,0,0.3,0.02]))
ax18.set_yticks([])
ax18.set_facecolor('black')
ax18.set_xticks([-100,0,100])
ax18.set_xticklabels([80,180,280])

plt.gcf().text(0.28,0.03, 'Mean Biomass Change Over 0-240m', fontsize=12)
plt.gcf().text(0.17,0.91, 'Kremer', fontsize=14, weight='bold')
plt.gcf().text(0.36,0.91, 'Anderson', fontsize=14, weight='bold')
plt.gcf().text(0.525,0.91, 'Anderson-Kremer', fontsize=14, weight='bold')
#plt.savefig('figures/FigureS6.pdf', bbox_inches='tight', transparent=True)

#%% Average change in biomass
# Figure 7

krem_diff = pft_sum_SQDUa-pft_sum_SQDU
and_diff = pft_sum_DQDUa-pft_sum_DQDU

diff_med_SQDU = np.nanmean(krem_diff*12.01/d, axis=2)
diff_med_DQDU = np.nanmean(and_diff*12.01/d, axis=2)

bins = np.linspace(0, 90, 7)

k_zmean = np.empty([6,6])
a_zmean = np.empty([6,6])

for i in range(6):
    for j in range(6):
        k_zmean[i,j] = np.nanmean(diff_med_SQDU[i , int(bins[j]):int(bins[j+1])])*12.01/d
        a_zmean[i,j] = np.nanmean(diff_med_DQDU[i , int(bins[j]):int(bins[j+1])])*12.01/d

fig, axs = plt.subplots(1,6, figsize=(12,4), sharex=True, sharey=True)
fig.subplots_adjust(hspace = .3, wspace=.15)

axs = axs.ravel()

left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
groups = ['Coccolithophores', 'Cyanobacteria', 'Diatoms', 'Diazotrophs','Dinoflagellates','Green Algae']
for i in range(6):
    k = k_zmean[i]
    a = a_zmean[i]
    index = ['[-90,-60]', '[-60,-30]', '[-30,0]', '[0,30]','[30,60]','[60,90]']
    ind = np.arange(6)
    width=0.4
    axs[i].barh(ind, k, width, color='#2A788EFF', label='Kremer')
    axs[i].barh(ind+width, a, width, color='#7AD151FF', label='Anderson')
    axs[i].set_title(groups[i], loc='center', fontsize=11)
    axs[0].set_ylabel('Latitude (Degrees)', fontsize=11)
    axs[0].set_yticks(ind + width / 2)
    axs[i].vlines(0,5.5,0, linestyle='--', color='grey', lw=1)
    axs[0].set_yticklabels(('[-90,-60]', '[-60,-30]', '[-30,0]', '[0,30]','[30,60]','[60,90]'))
plt.legend(bbox_to_anchor=(right+0.25, -.1))
plt.gcf().text(0.5 * (left + right), -0.06 * (bottom + top), 'Mean Change in Depth-Averaged Biomass \n(mg C $m^{-3}$)', 
               horizontalalignment='center',fontsize=11)

#plt.savefig('figures/Figure7.pdf', bbox_inches='tight', transparent=True)

