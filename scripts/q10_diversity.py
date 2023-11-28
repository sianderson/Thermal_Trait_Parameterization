#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Stephanie Anderson, Massachusetts Institute of Technology
Email: siander@mit.edu
    
This script contains calculations for phytoplankton community evenness and richness.

    INPUT: 
        Control simulations which has been archived here ############
        grid_igsm.nc: grid used in the model
    
    OUTPUT: 
        Figure 3 A&B. Phytoplankton community evenness and richness
        
%=========================================================================
"""

#%%
import os 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# depth calculations
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
     
#%% Depth integrated PP
# controls
# In all files, SQSU = Eppley, SQDU = Kremer, DQDU = Anderson
SQSUc_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_SQSU_control_avg.nc'
SQSUc = nc.Dataset(SQSUc_fname, 'r')

SQDUc_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_SQDU_control_avg.nc'
SQDUc = nc.Dataset(SQDUc_fname, 'r')

DQDUc_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_control_avg.nc'
DQDUc = nc.Dataset(DQDUc_fname, 'r')

# get primary production
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

#%% Average over top 260 m
# get primary production
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
    for i in ref[key]: # for each cell (tracer)
        cell = SQSUc.variables[i]
        cell2 = SQDUc.variables[i]
        cell3 = DQDUc.variables[i]
        cell_int = np.empty([8,90,144])   
        cell_sum = np.empty([9,90,144])
        cell_int2 = np.empty([8,90,144])   
        cell_sum2 = np.empty([9,90,144])
        cell_int3 = np.empty([8,90,144])   
        cell_sum3 = np.empty([9,90,144])
        for x in range(6): # 4 to 55m, 8 to 260m
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

#biomass_SQSU[biomass_SQSU < 0.001]=np.nan
#biomass_SQDU[biomass_SQDU < 0.001]=np.nan
#biomass_DQDU[biomass_DQDU < 0.001]=np.nan
#%% Evenness
s=5
import math

ref2 = ['TRAC43','TRAC42','TRAC41','TRAC40','TRAC39','TRAC38','TRAC37','TRAC36','TRAC35',
        'TRAC30','TRAC31','TRAC32','TRAC33','TRAC34',
        'TRAC25','TRAC26','TRAC27','TRAC28','TRAC29',
       'TRAC44','TRAC45','TRAC46','TRAC47','TRAC48','TRAC49','TRAC50','TRAC51',
       'TRAC23','TRAC24',
       'TRAC21','TRAC22']


# for each PFT
cell_biomass = np.empty([len(ref2),90,144])
cell_biomass2 = np.empty([len(ref2),90,144])
cell_biomass3 = np.empty([len(ref2),90,144])
num = 0
for i in ref2: 
    cell = SQSUc.variables[i]
    cell2 = SQDUc.variables[i]
    cell3 = DQDUc.variables[i]
    cell_int = np.empty([s,90,144])
    cell_int2 = np.empty([s,90,144])
    cell_int3 = np.empty([s,90,144])
    for x in range(s): # 4 to 55 m, 8 to 260 m
        cell_int[x] = depth_section[x] * cell[x]
        cell_int2[x] = depth_section[x] * cell2[x]
        cell_int3[x] = depth_section[x] * cell3[x]
    cell_biomass[num] = np.divide(np.sum(cell_int, axis=0), biomass_SQSU) 
    cell_biomass2[num] = np.divide(np.sum(cell_int2, axis=0), biomass_SQDU)
    cell_biomass3[num] = np.divide(np.sum(cell_int3, axis=0), biomass_DQDU)
    num = num+1

cell_biomass[cell_biomass < 0.001]=np.nan
cell_biomass2[cell_biomass2 < 0.001]=np.nan
cell_biomass3[cell_biomass3 < 0.001]=np.nan


H_SQSU = np.nansum(np.multiply(cell_biomass, np.log(cell_biomass)), axis=0)*-1
H_SQDU = np.nansum(np.multiply(cell_biomass2, np.log(cell_biomass2)), axis=0)*-1
H_DQDU = np.nansum(np.multiply(cell_biomass3, np.log(cell_biomass3)), axis=0)*-1

even_SQSU = H_SQSU/math.log(31)
even_SQDU = H_SQDU/math.log(31)
even_DQDU = H_DQDU/math.log(31)


#%% Species Richness

ref = {'diatoms':['TRAC43','TRAC42','TRAC41','TRAC40','TRAC39','TRAC38','TRAC37','TRAC36','TRAC35'],
       'diazotrophs':['TRAC30','TRAC31','TRAC32','TRAC33','TRAC34'],
       'coccolithophores':['TRAC25','TRAC26','TRAC27','TRAC28','TRAC29'],
       'dinos':['TRAC44','TRAC45','TRAC46','TRAC47','TRAC48','TRAC49','TRAC50','TRAC51'],
       'greens':['TRAC23','TRAC24'],
       'cyano':['TRAC21','TRAC22']}

pft_biomass_SQSU= np.empty([len(ref),90,144])
pft_biomass_SQDU = np.empty([len(ref),90,144])
pft_biomass_DQDU = np.empty([len(ref),90,144])
num1 = 0

for key in ref: # for each PFT
    cell_biomass = np.empty([len(ref[key]),90,144])
    cell_biomass2 = np.empty([len(ref[key]),90,144])
    cell_biomass3 = np.empty([len(ref[key]),90,144])
    num2 = 0
    for i in ref[key]: # for each cell (tracer)
        cell = SQSUc.variables[i]
        cell2 = SQDUc.variables[i]
        cell3 = DQDUc.variables[i]
        cell_int = np.empty([s,90,144])
        cell_int2 = np.empty([s,90,144])
        cell_int3 = np.empty([s,90,144])
        for x in range(s): # to 260 m
            cell_int[x] = depth_section[x] * cell[x]
            cell_int2[x] = depth_section[x] * cell2[x]
            cell_int3[x] = depth_section[x] * cell3[x]
        cell_biomass[num2] = np.divide(np.sum(cell_int, axis=0), biomass_SQSU)*100 # to 260 m
        cell_biomass2[num2] = np.divide(np.sum(cell_int2, axis=0), biomass_SQDU)*100
        cell_biomass3[num2] = np.divide(np.sum(cell_int3, axis=0), biomass_DQDU)*100
        num2 = num2+1
        cell_biomass_bi = np.where(cell_biomass[:] > 0.001, 1, 0)
        cell_biomass_bi2 = np.where(cell_biomass2[:] > 0.001, 1, 0)
        cell_biomass_bi3 = np.where(cell_biomass3[:] > 0.001, 1, 0) 
    pft_biomass_SQSU[num1] = np.sum(cell_biomass_bi, axis=0)
    pft_biomass_SQDU[num1] = np.sum(cell_biomass_bi2, axis=0)
    pft_biomass_DQDU[num1] = np.sum(cell_biomass_bi3, axis=0)
    num1 = num1+1
 
richness_SQSU = np.sum(pft_biomass_SQSU, axis=0)
richness_SQDU = np.sum(pft_biomass_SQDU, axis=0)
richness_DQDU = np.sum(pft_biomass_DQDU, axis=0)

#%% Richness and Evenness
land_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_control_avg.nc'
land = nc.Dataset(land_fname, 'r')
land = land.variables['PP']

# figure
import matplotlib
cmap2 = matplotlib.cm.get_cmap('viridis').copy()
cmap2.set_bad(color = 'k', alpha = 1)

even_SQSUm =np.ma.masked_where((biomass_SQSU < 0.001), even_SQSU)
even_SQDUm =np.ma.masked_where((biomass_SQDU < 0.001), even_SQDU)
even_DQDUm =np.ma.masked_where((biomass_DQDU < 0.001), even_DQDU)

richness_SQSUm =np.ma.masked_where((biomass_SQSU < 0.001), richness_SQSU)
richness_SQDUm =np.ma.masked_where((biomass_SQDU < 0.001), richness_SQDU)
richness_DQDUm =np.ma.masked_where((biomass_DQDU < 0.001), richness_DQDU)

gridspec = dict(hspace=-0.5,bottom=0, wspace=0.1, height_ratios=[1, 1, 1])
fig, [(ax1, ax4),(ax2, ax5), (ax3,ax6)] = plt.subplots(3, 2, figsize=(10, 9), 
                                                       sharey='row', sharex='col', gridspec_kw=gridspec)
im1 = ax1.imshow(even_SQSUm, extent=(0,360,-80,80),  vmin=0, vmax=1, cmap=cmap2)
ax1.invert_yaxis()
ax1.set_yticks([-50,0,50])
ax1.set_yticklabels([50,0,-50])
ax1.set_xticks([])
ax1.set_title('Eppley', loc='left')

im2= ax2.imshow(even_SQDUm, extent=(0,360,-80,80),  vmin=0, vmax=1, cmap=cmap2)
ax2.invert_yaxis()
ax2.set_yticks([-50,0,50])
ax2.set_yticklabels([50,0,-50])
ax2.set_xticks([])
ax2.set_title('Kremer', loc='left')

im3= ax3.imshow(even_DQDUm, extent=(0,360,-80,80),  vmin=0, vmax=1, cmap=cmap2)
ax3.invert_yaxis()
fig.colorbar(im3,label='Pielous Evenness', shrink=0.65, 
             orientation='horizontal',cax=fig.add_axes([0.16,0.06,0.3,0.02]))
ax3.set_yticks([-50,0,50])
ax3.set_xlabel('Longitude')
ax3.set_yticklabels([50,0,-50])
ax3.set_title('Anderson', loc='left')

im4 = ax4.imshow(richness_SQSUm, extent=(0,360,-80,80),  vmin=0, vmax=31, cmap=cmap2)
ax4.invert_yaxis()
ax4.set_yticks([-50,0,50])
ax4.set_yticklabels([50,0,-50])
ax4.set_xticks([])
ax4.set_title('Eppley', loc='left')

im5= ax5.imshow(richness_SQDUm, extent=(0,360,-80,80),  vmin=0, vmax=31, cmap=cmap2)
ax5.invert_yaxis()
ax5.set_yticks([-50,0,50])
ax5.set_yticklabels([50,0,-50])
ax5.set_xticks([])
ax5.set_title('Kremer', loc='left')

im6= ax6.imshow(richness_DQDUm, extent=(0,360,-80,80),  vmin=0, vmax=31, cmap=cmap2)
ax6.invert_yaxis()
fig.colorbar(im6,label='Species Richness', shrink=0.65, 
             orientation='horizontal',cax=fig.add_axes([0.56,0.06,0.3,0.02]))
ax6.set_yticks([-50,0,50])
ax6.set_xlabel('Longitude')
ax6.set_yticklabels([50,0,-50])
ax6.set_title('Anderson', loc='left')

plt.tight_layout()
fig.text(0.08,0.46, 'Latitude', va='center', rotation='vertical')
fig.text(0.09,0.77, '(a)',fontweight='bold')
fig.text(0.5,0.77, '(b)',fontweight='bold')
#plt.savefig('figures/Figure3AB.pdf', bbox_inches='tight', transparent=True) ## WRONG

#%% Richness/Evenness summary
# Figure 3A & 3B
rich_med_SQSU = np.nanmean(richness_SQSUm.filled(np.nan), axis=1)
rich_std_SQSU = np.nanstd(richness_SQSUm.filled(np.nan), axis=1) 
rich_med_SQDU = np.nanmean(richness_SQDUm.filled(np.nan), axis=1)
rich_std_SQDU = np.nanstd(richness_SQDUm.filled(np.nan), axis=1)
rich_med_DQDU = np.nanmean(richness_DQDUm.filled(np.nan), axis=1)
rich_std_DQDU = np.nanstd(richness_DQDUm.filled(np.nan), axis=1)

even_med_SQSU = np.nanmean(even_SQSUm.filled(np.nan), axis=1) 
even_std_SQSU = np.nanstd(even_SQSUm.filled(np.nan), axis=1) 
even_med_SQDU = np.nanmean(even_SQDUm.filled(np.nan), axis=1)
even_std_SQDU = np.nanstd(even_SQDUm.filled(np.nan), axis=1)
even_med_DQDU = np.nanmean(even_DQDUm.filled(np.nan), axis=1)
even_std_DQDU = np.nanstd(even_DQDUm.filled(np.nan), axis=1)

gridspec = dict(hspace=0.2,wspace=0.1,bottom=0)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3, 8), gridspec_kw=gridspec)
c2 = '#2A788EFF'
c3 = '#7AD151FF'

ax1.plot(rich_med_SQSU,range(-90,90,2), color='k')
ax1.fill_betweenx(range(-90,90,2),rich_med_SQSU-rich_std_SQSU, rich_med_SQSU+rich_std_SQSU, facecolor='k', alpha=0.1)
ax1.fill_betweenx(range(-90,90,2),rich_med_SQDU-rich_std_SQDU, rich_med_SQDU+rich_std_SQDU, facecolor=c2, alpha=0.1)
ax1.fill_betweenx(range(-90,90,2),rich_med_DQDU-rich_std_DQDU, rich_med_DQDU+rich_std_DQDU, facecolor=c3, alpha=0.1)
ax1.axhline(y=0, color='grey', linestyle='--')
ax1.plot(rich_med_SQSU,range(-90,90,2), color='k')
ax1.plot(rich_med_SQDU,range(-90,90,2), color=c2)
ax1.plot(rich_med_DQDU,range(-90,90,2), color=c3)
ax1.set_xlabel('Richness', fontsize=11)
ax1.set_ylabel('Latitude', fontsize=11)
ax1.set_xlim(7,25)

ax2.axhline(y=0, color='grey', linestyle='--')
line1, = ax2.plot(even_med_SQSU,range(-90,90,2), color='k',label='Eppley')
line2, = ax2.plot(even_med_SQDU,range(-90,90,2), color=c2,label='Kremer')
line3, = ax2.plot(even_med_DQDU,range(-90,90,2), color=c3,label='Anderson')
ax2.set_ylabel('Latitude', fontsize=11)
ax2.set_xlim(0.2,1)
ax2.fill_betweenx(range(-90,90,2),even_med_SQSU-even_std_SQSU, even_med_SQSU+even_std_SQSU, facecolor='k', alpha=0.1)
ax2.fill_betweenx(range(-90,90,2),even_med_SQDU-even_std_SQDU, even_med_SQDU+even_std_SQDU, facecolor=c2, alpha=0.1)
ax2.fill_betweenx(range(-90,90,2),even_med_DQDU-even_std_DQDU, even_med_DQDU+even_std_DQDU, facecolor=c3, alpha=0.1)
leg = ax2.legend(handles = [line1,line2,line3],bbox_to_anchor=(0.3, 0.08), loc="lower center", frameon=False)
for line in leg.get_lines():
    line.set_linewidth(3.0)
ax2.set_xlabel('Evenness', fontsize=11)

#plt.savefig("figures/Figure3AB.pdf", bbox_inches='tight', transparent=True)

