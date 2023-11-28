#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Stephanie Anderson, Massachusetts Institute of Technology
Email: siander@mit.edu
    
This script examines community change using Bray-Curtis dissimilarity.

    INPUT: 
        Control simulations which has been archived here ############
        grid_igsm.nc: grid used in the model
        Requires PFT biomass from control simulations in q10_bray_curtis.py
    
    OUTPUT: 
       Figure S9. Primary production temp only vs. all physics changing
       Figure S11. Bray-Curtis dissimilarity
     
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
# In all files, SQSU = Eppley, SQDU = Kremer, DQDU = Anderson

# Controls 
SQSUc_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_SQSU_control_avg.nc'
SQSUc = nc.Dataset(SQSUc_fname, 'r')

SQDUc_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_SQDU_control_avg.nc'
SQDUc = nc.Dataset(SQDUc_fname, 'r')

DQDUc_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_control_avg.nc'
DQDUc = nc.Dataset(DQDUc_fname, 'r')

# Climate Change (all physics changing)
SQSU_2080all_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_SQSU_2080avg.nc'
SQSUall = nc.Dataset(SQSU_2080all_fname, 'r')

SQDU_2080all_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_SQDU_2080avg.nc'
SQDUall = nc.Dataset(SQDU_2080all_fname, 'r') 

DQDU_2080all_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_2080avg.nc'
DQDUall = nc.Dataset(DQDU_2080all_fname, 'r')

# Climate Change (warming only)
SQSU_2080temp_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_SQSU_2080temp_avg.nc'
SQSUtemp = nc.Dataset(SQSU_2080temp_fname, 'r')

SQDU_2080temp_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_SQDU_2080temp_avg.nc'
SQDUtemp = nc.Dataset(SQDU_2080temp_fname, 'r')

DQDU_2080temp_fname = '/Volumes/SABackup/MIT_q10/Ensembles/run33_5_ice_2080temp_avg.nc'
DQDUtemp = nc.Dataset(DQDU_2080temp_fname, 'r')

# get primary production
PP_SQSUc = SQSUc.variables['PP'] 
PP_SQDUc = SQDUc.variables['PP']
PP_DQDUc = DQDUc.variables['PP']

PP_SQSUall = SQSUall.variables['PP']
PP_SQDUall = SQDUall.variables['PP']
PP_DQDUall = DQDUall.variables['PP']

PP_SQSUtemp = SQSUtemp.variables['PP']
PP_SQDUtemp = SQDUtemp.variables['PP']
PP_DQDUtemp = DQDUtemp.variables['PP']

PP_int_SQSUc = np.empty([22,90,144])
PP_int_sum_SQSUc = np.empty([22,90,144])
PP_int_SQDUc = np.empty([22,90,144])
PP_int_sum_SQDUc = np.empty([22,90,144])
PP_int_DQDUc = np.empty([22,90,144])
PP_int_sum_DQDUc = np.empty([22,90,144])

PP_int_SQSUall = np.empty([22,90,144])
PP_int_sum_SQSUall = np.empty([22,90,144])
PP_int_SQDUall = np.empty([22,90,144])
PP_int_sum_SQDUall = np.empty([22,90,144])
PP_int_DQDUall = np.empty([22,90,144])
PP_int_sum_DQDUall = np.empty([22,90,144])

PP_int_SQSUt = np.empty([22,90,144])
PP_int_sum_SQSUt = np.empty([22,90,144])
PP_int_SQDUt = np.empty([22,90,144])
PP_int_sum_SQDUt = np.empty([22,90,144])
PP_int_DQDUt = np.empty([22,90,144])
PP_int_sum_DQDUt = np.empty([22,90,144])

for x in range(22):
    PP_int_SQSUc[x] = depth_section[x] * PP_SQSUc[x]
    PP_int_SQDUc[x] = depth_section[x] * PP_SQDUc[x]
    PP_int_DQDUc[x] = depth_section[x] * PP_DQDUc[x]
    PP_int_SQSUall[x] = depth_section[x] * PP_SQSUall[x]
    PP_int_SQDUall[x] = depth_section[x] * PP_SQDUall[x]
    PP_int_DQDUall[x] = depth_section[x] * PP_DQDUall[x]
    PP_int_SQSUt[x] = depth_section[x] * PP_SQSUtemp[x]
    PP_int_SQDUt[x] = depth_section[x] * PP_SQDUtemp[x]
    PP_int_DQDUt[x] = depth_section[x] * PP_DQDUtemp[x]
    if x<1:
        PP_int_sum_SQSUc[x]=PP_int_SQSUc[x]
        PP_int_sum_SQDUc[x]=PP_int_SQDUc[x]
        PP_int_sum_DQDUc[x]=PP_int_DQDUc[x]
        PP_int_sum_SQSUall[x]=PP_int_SQSUall[x]
        PP_int_sum_SQDUall[x]=PP_int_SQDUall[x]
        PP_int_sum_DQDUall[x]=PP_int_DQDUall[x]
        PP_int_sum_SQSUt[x]=PP_int_SQSUt[x]
        PP_int_sum_SQDUt[x]=PP_int_SQDUt[x]
        PP_int_sum_DQDUt[x]=PP_int_DQDUt[x]
    if x>=1:
        PP_int_sum_SQSUc[x]=PP_int_sum_SQSUc[x-1]+PP_int_SQSUc[x]
        PP_int_sum_SQDUc[x]=PP_int_sum_SQDUc[x-1]+PP_int_SQDUc[x]
        PP_int_sum_DQDUc[x]=PP_int_sum_DQDUc[x-1]+PP_int_DQDUc[x]
        PP_int_sum_SQSUall[x]=PP_int_sum_SQSUall[x-1]+PP_int_SQSUall[x]
        PP_int_sum_SQDUall[x]=PP_int_sum_SQDUall[x-1]+PP_int_SQDUall[x]
        PP_int_sum_DQDUall[x]=PP_int_sum_DQDUall[x-1]+PP_int_DQDUall[x]
        PP_int_sum_SQSUt[x]=PP_int_sum_SQSUt[x-1]+PP_int_SQSUt[x]
        PP_int_sum_SQDUt[x]=PP_int_sum_SQDUt[x-1]+PP_int_SQDUt[x]
        PP_int_sum_DQDUt[x]=PP_int_sum_DQDUt[x-1]+PP_int_DQDUt[x]
  

#%% Difference in PP with Warming only vs. all physics changing (Figure S9)
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors

c = 0.01201*3.154e+7 # (adjusting units to moles per year)
s = 6 # Depth (240 m)
vmin = -40
vmax = 40
PP_SQSU_ctemp = np.subtract(PP_int_sum_SQSUt*c, PP_int_sum_SQSUc*c)
PP_SQDU_ctemp = np.subtract(PP_int_sum_SQDUt*c, PP_int_sum_SQDUc*c)
PP_DQDU_ctemp = np.subtract(PP_int_sum_DQDUt*c, PP_int_sum_DQDUc*c)

PP_SQSU_call= np.subtract(PP_int_sum_SQSUall*c, PP_int_sum_SQSUc*c)
PP_SQDU_call = np.subtract(PP_int_sum_SQDUall*c, PP_int_sum_SQDUc*c)
PP_DQDU_call = np.subtract(PP_int_sum_DQDUall*c, PP_int_sum_DQDUc*c)

pp_med_SQSU_ctemp = np.nanmean(PP_SQSU_ctemp[s], axis=1)
pp_std_SQSU_ctemp = np.nanstd(PP_SQSU_ctemp[s], axis=1) 
pp_med_SQDU_ctemp = np.nanmean(PP_SQDU_ctemp[s], axis=1)
pp_std_SQDU_ctemp = np.nanstd(PP_SQDU_ctemp[s], axis=1)
pp_med_DQDU_ctemp = np.nanmean(PP_DQDU_ctemp[s], axis=1)
pp_std_DQDU_ctemp = np.nanstd(PP_DQDU_ctemp[s], axis=1)

pp_med_SQSU_call = np.nanmean(PP_SQSU_call[s], axis=1)
pp_std_SQSU_call = np.nanstd(PP_SQSU_call[s], axis=1) 
pp_med_SQDU_call = np.nanmean(PP_SQDU_call[s], axis=1)
pp_std_SQDU_call = np.nanstd(PP_SQDU_call[s], axis=1)
pp_med_DQDU_call = np.nanmean(PP_DQDU_call[s], axis=1)
pp_std_DQDU_call = np.nanstd(PP_DQDU_call[s], axis=1)



cmap2 = matplotlib.cm.get_cmap('RdBu_r').copy()
cmap2.set_bad(color = 'white', alpha = 1)

grid = dict(hspace=0.3,bottom=0.1,height_ratios=[1,1],wspace=0.04, width_ratios=[1,1,1,0.45])
fig, [(ax1, ax2,ax3, ax7), (ax4,ax5,ax6, ax8)] = plt.subplots(2, 4, figsize=(12, 4), 
                                                              sharex='col',
                                                              sharey='row',
                                                              gridspec_kw=grid,
                                                              subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180, globe=None)})

im1 = ax1.imshow(PP_SQSU_ctemp[21], extent=(0,360,-90,90), cmap=cmap2, vmin=vmin, vmax=vmax, 
                 origin='lower', transform=ccrs.PlateCarree())
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax1.set_title('Temperature Only', loc='left')
ax1.set_yticks([-50,0,50])

im2= ax2.imshow(PP_SQDU_ctemp[21], extent=(0,360,-90,90), cmap=cmap2, vmin=vmin, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

im3= ax3.imshow(PP_DQDU_ctemp[21], extent=(0,360,-90,90), cmap=cmap2, vmin=vmin, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax3.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

im4 = ax4.imshow(PP_SQSU_call[21], extent=(0,360,-90,90), cmap=cmap2, vmin=vmin, vmax=vmax,
                 origin='lower', transform=ccrs.PlateCarree())
ax4.set_yticks([-50,0,50])
ax4.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax4.set_title('All Physics', loc='left')
ax4.set_xticks([-180,-80,20,120])
ax4.set_xticklabels(['0','100','200','300'])

im5= ax5.imshow(PP_SQDU_call[21], extent=(0,360,-90,90), cmap=cmap2, vmin=vmin, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax5.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax5.set_xticks([-180,-80,20,120])
ax5.set_xticklabels(['0','100','200','300'])

im6= ax6.imshow(PP_DQDU_call[21], extent=(0,360,-90,90), cmap=cmap2, vmin=vmin, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
fig.colorbar(im3,label='PP Change (g C $m^{-2}$ $y^{-1}$)', shrink=0.65, 
             orientation='horizontal', norm=colors.CenteredNorm(),cax=fig.add_axes([0.33,-0.08,0.3,0.03]))
ax6.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax6.set_xticks([-180,-80,20,120])
ax6.set_xticklabels(['0','100','200','300'])
ax5.set_xlabel('Longitude')

c2 = '#2A788EFF'
c3 = '#7AD151FF'

r= 2
ax7.plot(pp_med_SQSU_ctemp*r,range(-90,90,2), color='k')
ax7.axvline(0, color='gray', linestyle='--')
ax7.fill_betweenx(range(-90,90,2),(pp_med_SQSU_ctemp-pp_std_SQSU_ctemp)*r, (pp_med_SQSU_ctemp+pp_std_SQSU_ctemp)*r, facecolor='k', alpha=0.1)
ax7.fill_betweenx(range(-90,90,2),(pp_med_SQDU_ctemp-pp_std_SQDU_ctemp)*r, (pp_med_SQDU_ctemp+pp_std_SQDU_ctemp)*r, facecolor=c2, alpha=0.1)
ax7.fill_betweenx(range(-90,90,2),(pp_med_DQDU_ctemp-pp_std_DQDU_ctemp)*r, (pp_med_DQDU_ctemp+pp_std_DQDU_ctemp)*r, facecolor=c3, alpha=0.1)
ax7.plot(pp_med_SQSU_ctemp*r,range(-90,90,2), color='k')
ax7.plot(pp_med_SQDU_ctemp*r,range(-90,90,2), color=c2)
ax7.plot(pp_med_DQDU_ctemp*r,range(-90,90,2), color=c3)
ax7.set_xlim([-80, 80])

line1, = ax8.plot(pp_med_SQSU_call*r,range(-90,90,2), color='k',label='Eppley')
ax8.axvline(0, color='gray', linestyle='--')
line2, = ax8.plot(pp_med_SQDU_call*r,range(-90,90,2), color=c2,label='Kremer')
line3, = ax8.plot(pp_med_DQDU_call*r,range(-90,90,2), color=c3,label='Anderson')
ax8.fill_betweenx(range(-90,90,2),(pp_med_SQSU_call-pp_std_SQSU_call)*r, (pp_med_SQSU_call+pp_std_SQSU_call)*r, facecolor='k', alpha=0.1)
ax8.fill_betweenx(range(-90,90,2),(pp_med_SQDU_call-pp_std_SQDU_call)*r, (pp_med_SQDU_call+pp_std_SQDU_call)*r, facecolor=c2, alpha=0.1)
ax8.fill_betweenx(range(-90,90,2),(pp_med_DQDU_call-pp_std_DQDU_call)*r, (pp_med_DQDU_call+pp_std_DQDU_call)*r, facecolor=c3, alpha=0.1)
leg = ax8.legend(handles = [line1,line2,line3],bbox_to_anchor=(0.5, -0.9), loc="lower center", frameon=False)
for line in leg.get_lines():
    line.set_linewidth(3.0)
ax8.set_xlabel('Zonal Mean', fontsize=11)
ax8.set_xticks([-60,0,60])
ax8.set_xticklabels([-30,0,30])

plt.gcf().text(0.20,0.96, 'Eppley', fontsize=14, weight='bold')
plt.gcf().text(0.44,0.96, 'Kremer', fontsize=14, weight='bold')
plt.gcf().text(0.66,0.96, 'Anderson', fontsize=14, weight='bold')
plt.gcf().text(0.07,0.45, 'Latitude', fontsize=11, rotation=90)
#plt.savefig('figures/FigureS9.pdf', bbox_inches='tight', transparent=True)


#%% Average over top 260 m
ref = {'coccolithophores':['TRAC25','TRAC26','TRAC27','TRAC28','TRAC29'],
       'cyano':['TRAC21','TRAC22'],
       'diatoms':['TRAC43','TRAC42','TRAC41','TRAC40','TRAC39','TRAC38','TRAC37','TRAC36','TRAC35'],
       'diazotrophs':['TRAC30','TRAC31','TRAC32','TRAC33','TRAC34'],
       'dinos':['TRAC44','TRAC45','TRAC46','TRAC47','TRAC48','TRAC49','TRAC50','TRAC51'],
       'greens':['TRAC23','TRAC24'],
       }

pft_sum_SQDUt = np.empty([len(ref),90,144])
pft_sum_DQDUt = np.empty([len(ref),90,144])
num1 = 0

for key in ref: # for each PFT
    all_cell_sum2 = np.empty([len(ref[key]),90,144])
    all_cell_sum3 = np.empty([len(ref[key]),90,144])
    num2 = 0
    for i in ref[key]: # for each phenotype (tracer)
        cell2 = SQDUtemp.variables[i]
        cell3 = DQDUtemp.variables[i]
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
    pft_sum_SQDUt[num1] = np.sum(all_cell_sum2, axis=0)
    pft_sum_DQDUt[num1] = np.sum(all_cell_sum3, axis=0)
    num1 = num1+1

biomass_DQDUt = np.nansum(pft_sum_DQDUt, axis=0)
biomass_SQDUt = np.nansum(pft_sum_SQDUt, axis=0)

biomass_SQDUt[biomass_SQDUt < 0.001]=np.nan
biomass_DQDUt[biomass_DQDUt < 0.001]=np.nan

pft_biomass_SQDUt = np.empty([len(ref),90,144])
pft_biomass_DQDUt = np.empty([len(ref),90,144])
for i in range(6):
    pft_biomass_SQDUt[i] = np.divide(pft_sum_SQDUt[i],biomass_SQDUt)
    pft_biomass_DQDUt[i] = np.divide(pft_sum_DQDUt[i],biomass_DQDUt)
    
#%% Change in community structure (PFT)
# Requires controls from q10_bray_curtis.py
# (pft_biomass_DQDU & pft_biomass_SQDU)
     
c2 = np.empty([6])
d2 = np.empty([6])
c3 = np.empty([6])
d3 = np.empty([6])

c4 = np.empty([6])
d4 = np.empty([6])
c5 = np.empty([6])
d5 = np.empty([6])

SQSU_brayt = np.empty([90,144])
SQDU_brayt = np.empty([90,144])
DQDU_brayt = np.empty([90,144])
bray_control = np.empty([90,144])
bray_t = np.empty([90,144])

for j in range(90):
    for k in range(144):
        for i in range(6):
            c2[i]=abs(pft_biomass_SQDUt[i,j,k] - pft_biomass_SQDU[i,j,k])
            d2[i]=pft_biomass_SQDUt[i,j,k] + pft_biomass_SQDU[i,j,k]
            c3[i]=abs(pft_biomass_DQDUt[i,j,k] - pft_biomass_DQDU[i,j,k])
            d3[i]=pft_biomass_DQDUt[i,j,k] + pft_biomass_DQDU[i,j,k]
            
            c4[i]=abs(pft_biomass_DQDU[i,j,k] - pft_biomass_SQDU[i,j,k])
            d4[i]=pft_biomass_DQDU[i,j,k] + pft_biomass_SQDU[i,j,k]
            c5[i]=abs(pft_biomass_DQDUt[i,j,k] - pft_biomass_SQDUt[i,j,k])
            d5[i]=pft_biomass_DQDUt[i,j,k] + pft_biomass_SQDUt[i,j,k]
            
      
        SQDU_brayt[j,k]=np.divide(np.nansum(c2),np.nansum(d2),
                                    out=np.zeros_like(np.nansum(c2)), where=np.nansum(d2)!=0)
        DQDU_brayt[j,k]=np.divide(np.nansum(c3),np.nansum(d3),
                                     out=np.zeros_like(np.nansum(c3)), where=np.nansum(d3)!=0)
        
        bray_control[j,k]=np.divide(np.nansum(c4),np.nansum(d4),
                                   out=np.zeros_like(np.nansum(c4)), where=np.nansum(d4)!=0)
        bray_t[j,k]=np.divide(np.nansum(c5),np.nansum(d5),
                                     out=np.zeros_like(np.nansum(c5)), where=np.nansum(d5)!=0)


#%% Figure S11 (Bray-Curtis: Warming ONLY)
import matplotlib
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

cmap2 = matplotlib.cm.get_cmap('cividis').copy()
cmap2.set_bad(color = 'k', alpha = 1)
vmax=0.4

gridspec = dict(height_ratios=[1,1],width_ratios=[1,1], hspace=0.2, wspace=0)
fig, [(ax5,ax7),(ax6, ax8)] = plt.subplots(2,2, figsize=(8,4), sharey='row', gridspec_kw=gridspec)

ax5 = plt.subplot2grid((2, 2), (0, 0), projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im5= ax5.imshow(SQDU_brayt, extent=(0,360,-90,90), cmap=cmap2, vmin=0, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax5.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax5.set_yticks([-50,0,50])
ax5.set_title('Kremer',loc='left')

ax6 = plt.subplot2grid((2, 2), (1, 0), projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im6= ax6.imshow(DQDU_brayt, extent=(0,360,-90,90), cmap=cmap2, vmin=0, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax6.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax6.set_title('Anderson',loc='left')
ax6.set_xticks([-100,0,100])
ax6.set_xticklabels([80,180,280])
ax6.set_yticks([-50,0,50])

ax7 = plt.subplot2grid((2, 2), (0, 1), projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im7= ax7.imshow(bray_control, extent=(0,360,-90,90), cmap=cmap2, vmin=0, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax7.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax7.set_title('Control',loc='left')


ax8 = plt.subplot2grid((2, 2), (1, 1), projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im8= ax8.imshow(bray_t, extent=(0,360,-90,90), cmap=cmap2, vmin=0, vmax=vmax,
                origin='lower', transform=ccrs.PlateCarree())
ax8.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im8,label='Bray-Curtis Dissimilarity', shrink=0.65, 
             orientation='horizontal', norm=colors.CenteredNorm(),cax=fig.add_axes([0.37,-0.04,0.3,0.03]))
ax8.set_title('Warming Only',loc='left')
ax8.set_xticks([-100,0,100])
ax8.set_xticklabels([80,180,280])

fig.text(0.07,0.5, 'Latitude', va='center', rotation='vertical', size=12)
fig.text(0.45,0.04, 'Longitude', va='center', size=12)
plt.gcf().text(0.15,0.95, '1860 vs 2100 (Warming Only)', fontsize=12, weight='bold')
plt.gcf().text(0.535,0.95, 'Kremer vs Anderson', fontsize=12, weight='bold')

#plt.savefig('figures/FigureS11.pdf', bbox_inches='tight', transparent=True)
