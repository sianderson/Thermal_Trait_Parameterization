#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Stephanie Anderson, Massachusetts Institute of Technology
Email: siander@mit.edu
    
This script examines primary production (PP) and biomass in the model runs. 

    INPUT: 
        Control simulations which have been archived here: https://doi.org/10.7910/DVN/6TLL8Z
        grid_igsm.nc: grid used in the model
    
    OUTPUT: 
        Figure 2. Biomass and PP   
        Figure 4. Mean biomass
        Figure S5. Change in PP and biomass
        
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
depth2 = np.empty([23])
depth_section = np.empty([22])
for x in range(22):
    if x==0:
        depth2[x] = 2*depth[x]
        depth_section[x] = depth2[x]
    if x>0:
        depth2[x] = 2*depth[x]-depth2[x-1]
        depth_section[x] = depth2[x] - depth2[x-1]
     
#%% Depth integrated PP
# Controls
# In all files, SQSU = Eppley, SQDU = Kremer, DQDU = Anderson
SQSUc_fname = 'Eppley_control.nc'
SQSUc = nc.Dataset(SQSUc_fname, 'r')

SQDUc_fname = 'Kremer_control.nc'
SQDUc = nc.Dataset(SQDUc_fname, 'r')

DQDUc_fname = 'Anderson_control.nc'
DQDUc = nc.Dataset(DQDUc_fname, 'r')

# get primary production
PP_SQSUc = SQSUc.variables['PP'] # units  =  mmol C /m3 / s
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
        

#%% difference in PP over 240 m
import matplotlib
import matplotlib.colors as colors

s=6
PP_SQSUc_SQDUc = np.divide(np.subtract(PP_int_sum_SQDUc[s], PP_int_sum_SQSUc[s]),
                               PP_int_sum_SQSUc[s])*100

PP_SQSUc_DQDUc = np.divide(np.subtract(PP_int_sum_DQDUc[s], PP_int_sum_SQSUc[s]),
                               PP_int_sum_SQSUc[s])*100

PP_SQDUc_DQDUc = np.divide(np.subtract(PP_int_sum_DQDUc[s], PP_int_sum_SQDUc[s]),
                               PP_int_sum_SQDUc[s])*100

# mean difference
np.nanmean([PP_SQSUc_SQDUc[0:15,],PP_SQSUc_SQDUc[75:90,]])
np.nanmean([PP_SQSUc_DQDUc[0:15,],PP_SQSUc_DQDUc[75:90,]])
np.nanmean([np.nanmean([PP_SQSUc_SQDUc[0:15,],PP_SQSUc_SQDUc[75:90,]]),
            np.nanmean([PP_SQSUc_DQDUc[0:15,],PP_SQSUc_DQDUc[75:90,]])])

#%% Average over top 240
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
        cell = SQSUc.variables[i] # mmol C /m3
        cell2 = SQDUc.variables[i]
        cell3 = DQDUc.variables[i]
        cell_int = np.empty([8,90,144])   
        cell_sum = np.empty([9,90,144])
        cell_int2 = np.empty([8,90,144])   
        cell_sum2 = np.empty([9,90,144])
        cell_int3 = np.empty([8,90,144])   
        cell_sum3 = np.empty([9,90,144])
        for x in range(6): # 6 = 240 m
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


#%% PFT Biomass (Figure 4)
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

d=240 # depth
vmax=10
pft_med_SQSU = np.nanmean(pft_sum_SQSU*12.01/d, axis=2)
pft_std_SQSU = np.nanstd(pft_sum_SQSU*12.01/d, axis=2) 
pft_med_SQDU = np.nanmean(pft_sum_SQDU*12.01/d, axis=2)
pft_std_SQDU = np.nanstd(pft_sum_SQDU*12.01/d, axis=2)
pft_med_DQDU = np.nanmean(pft_sum_DQDU*12.01/d, axis=2)
pft_std_DQDU = np.nanstd(pft_sum_DQDU*12.01/d, axis=2)

cmap = matplotlib.cm.get_cmap('viridis')
cmap.set_bad(color = '#440154', alpha = 1)

gridspec = dict(hspace=0.1,bottom=0.05, width_ratios=[1,1,1,0.5], wspace=0.1)
fig, ([ax1, ax7, ax13, ax19], [ax2, ax8,ax14, ax20], [ax3, ax9, ax15, ax21],
      [ax4,ax10,ax16, ax22], [ax5,ax11, ax17, ax23], [ax6,ax12,ax18, ax24]) = plt.subplots(6,4, figsize=(11, 10), sharey='col', sharex='col',gridspec_kw=gridspec)
                                                                                         
ax1 = plt.subplot(641,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im1 = ax1.imshow((pft_sum_SQSU[0]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                 origin='lower', transform=ccrs.PlateCarree())
ax1.set_yticks([-50,0,50])
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax1.set_facecolor('black')
ax1.set_title('Coccolithophores', loc='left')

ax2 = plt.subplot(645,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im2= ax2.imshow((pft_sum_SQSU[1]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                origin='lower', transform=ccrs.PlateCarree())
ax2.set_yticks([-50,0,50])
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax2.set_facecolor('black')
ax2.set_title('Cyanobacteria', loc='left')

ax3 = plt.subplot(649,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im3= ax3.imshow((pft_sum_SQSU[2]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                origin='lower', transform=ccrs.PlateCarree())
ax3.set_yticks([-50,0,50])
ax3.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax3.set_facecolor('black')
ax3.set_title('Diatoms', loc='left')

ax4 = plt.subplot(6,4,13,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im4= ax4.imshow((pft_sum_SQSU[3]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
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
im5= ax5.imshow((pft_sum_SQSU[4]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                origin='lower', transform=ccrs.PlateCarree())
ax5.set_yticks([-50,0,50])
ax5.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax5.set_facecolor('black')
ax5.set_title('Dinoflagellates', loc='left')

ax6 = plt.subplot(6,4,21,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im6= ax6.imshow((pft_sum_SQSU[5]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                origin='lower', transform=ccrs.PlateCarree())
ax6.set_yticks([-50,0,50])
ax6.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax6.set_facecolor('black')
ax6.set_title('Green Algae', loc='left')
ax6.set_xticks([-100,0,100])
ax6.set_xticklabels([80,180,280])

ax7 = plt.subplot(6,4,2,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im7 = ax7.imshow((pft_sum_SQDU[0]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                 origin='lower', transform=ccrs.PlateCarree())#norm=matplotlib.colors.LogNorm())
ax7.set_yticks([])
ax7.set_facecolor('black')
ax7.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))

ax8 = plt.subplot(6,4,6,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im8= ax8.imshow((pft_sum_SQDU[1]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                origin='lower', transform=ccrs.PlateCarree())
ax8.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax8.set_yticks([])
ax8.set_facecolor('black')

ax9 = plt.subplot(6,4,10,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im9= ax9.imshow((pft_sum_SQDU[2]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                origin='lower', transform=ccrs.PlateCarree())
ax9.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax9.set_yticks([])
ax9.set_facecolor('black')

ax10 = plt.subplot(6,4,14,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im10= ax10.imshow((pft_sum_SQDU[3]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                  origin='lower', transform=ccrs.PlateCarree())
ax10.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax10.set_facecolor('black')
ax10.set_yticks([])

ax11 = plt.subplot(6,4,18,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im11= ax11.imshow((pft_sum_SQDU[4]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                  origin='lower', transform=ccrs.PlateCarree())
ax11.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax11.set_yticks([])
ax11.set_facecolor('black')

ax12 = plt.subplot(6,4,22,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im12= ax12.imshow((pft_sum_SQDU[5]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                  origin='lower', transform=ccrs.PlateCarree())
ax12.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax12.set_yticks([])
ax12.set_xlabel('Longitude', fontsize=12)
ax12.set_xticks([-100,0,100])
ax12.set_xticklabels([80,180,280])

ax13 = plt.subplot(6,4,3,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im13 = ax13.imshow((pft_sum_DQDU[0]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                   origin='lower', transform=ccrs.PlateCarree())
ax13.set_yticks([])
ax13.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax13.set_facecolor('black')

ax14 = plt.subplot(6,4,7,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im14= ax14.imshow((pft_sum_DQDU[1]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                  origin='lower', transform=ccrs.PlateCarree())
ax14.set_yticks([])
ax14.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax14.set_yticks([])
ax14.set_facecolor('black')

ax15 = plt.subplot(6,4,11,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im15= ax15.imshow((pft_sum_DQDU[2]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                  origin='lower', transform=ccrs.PlateCarree())
ax15.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax15.set_yticks([])
ax15.set_facecolor('black')

ax16 = plt.subplot(6,4,15,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im16= ax16.imshow((pft_sum_DQDU[3]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                  origin='lower', transform=ccrs.PlateCarree())
ax16.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax16.set_facecolor('black')
ax16.set_yticks([])

ax17 = plt.subplot(6,4,19,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im17= ax17.imshow((pft_sum_DQDU[4]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                  origin='lower', transform=ccrs.PlateCarree())
ax17.set_yticks([])
ax17.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax17.set_facecolor('black')

ax18 = plt.subplot(6,4,23,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im18= ax18.imshow((pft_sum_DQDU[5]*12.01/d), extent=(0,360,-90,90), norm=LogNorm(vmin=0.1, vmax=vmax), cmap=cmap,
                  origin='lower', transform=ccrs.PlateCarree())
ax18.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im18,label='mg C $m^{-3}$', orientation="horizontal",  cax=fig.add_axes([0.3,-0.04,0.3,0.02]))
ax18.set_yticks([])
ax18.set_facecolor('black')
ax18.set_xticks([-100,0,100])
ax18.set_xticklabels([80,180,280])

c2 = '#2A788EFF'
c3 = '#7AD151FF'
ax19 = plt.subplot(6,4,4)
ax19.plot(pft_med_SQSU[0],range(-90,90,2), color='k')
ax19.plot(pft_med_SQDU[0],range(-90,90,2), color=c2)
ax19.plot(pft_med_DQDU[0],range(-90,90,2), color=c3)
ax19.set_yticks([])
ax19.fill_betweenx(range(-90,90,2),pft_med_SQSU[0]-pft_std_SQSU[0], pft_med_SQSU[0]+pft_std_SQSU[0], facecolor='k', alpha=0.2)
ax19.fill_betweenx(range(-90,90,2),pft_med_SQDU[0]-pft_std_SQDU[0], pft_med_SQDU[0]+pft_std_SQDU[0], facecolor=c2, alpha=0.3)
ax19.fill_betweenx(range(-90,90,2),pft_med_DQDU[0]-pft_std_DQDU[0], pft_med_DQDU[0]+pft_std_DQDU[0], facecolor=c3, alpha=0.3)

ax20 = plt.subplot(6,4,8)
ax20.plot(pft_med_SQSU[1],range(-90,90,2), color='k')
ax20.plot(pft_med_SQDU[1],range(-90,90,2), color=c2)
ax20.plot(pft_med_DQDU[1],range(-90,90,2), color=c3)
ax20.fill_betweenx(range(-90,90,2),pft_med_SQSU[1]-pft_std_SQSU[1], pft_med_SQSU[1]+pft_std_SQSU[1], facecolor='k', alpha=0.2)
ax20.fill_betweenx(range(-90,90,2),pft_med_SQDU[1]-pft_std_SQDU[1], pft_med_SQDU[1]+pft_std_SQDU[1], facecolor=c2, alpha=0.3)
ax20.fill_betweenx(range(-90,90,2),pft_med_DQDU[1]-pft_std_DQDU[1], pft_med_DQDU[1]+pft_std_DQDU[1], facecolor=c3, alpha=0.3)

ax21 = plt.subplot(6,4,12)
ax21.plot(pft_med_SQSU[2],range(-90,90,2), color='k')
ax21.plot(pft_med_SQDU[2],range(-90,90,2), color='#39568CFF')
ax21.plot(pft_med_DQDU[2],range(-90,90,2), color='#3CBB75FF')
ax21.fill_betweenx(range(-90,90,2),pft_med_SQSU[2]-pft_std_SQSU[2], pft_med_SQSU[2]+pft_std_SQSU[2], facecolor='k', alpha=0.2)
ax21.fill_betweenx(range(-90,90,2),pft_med_SQDU[2]-pft_std_SQDU[2], pft_med_SQDU[2]+pft_std_SQDU[2], facecolor=c2, alpha=0.3)
ax21.fill_betweenx(range(-90,90,2),pft_med_DQDU[2]-pft_std_DQDU[2], pft_med_DQDU[2]+pft_std_DQDU[2], facecolor=c3, alpha=0.3)

ax22 = plt.subplot(6,4,16)
ax22.plot(pft_med_SQSU[3],range(-90,90,2), color='k')
ax22.plot(pft_med_SQDU[3],range(-90,90,2), color=c2)
ax22.plot(pft_med_DQDU[3],range(-90,90,2), color=c3)
ax22.fill_betweenx(range(-90,90,2),pft_med_SQSU[3]-pft_std_SQSU[3], pft_med_SQSU[3]+pft_std_SQSU[3], facecolor='k', alpha=0.2)
ax22.fill_betweenx(range(-90,90,2),pft_med_SQDU[3]-pft_std_SQDU[3], pft_med_SQDU[3]+pft_std_SQDU[3], facecolor=c2, alpha=0.3)
ax22.fill_betweenx(range(-90,90,2),pft_med_DQDU[3]-pft_std_DQDU[3], pft_med_DQDU[3]+pft_std_DQDU[3], facecolor=c3, alpha=0.3)

ax23 = plt.subplot(6,4,20)
ax23.plot(pft_med_SQSU[4],range(-90,90,2), color='k')
ax23.plot(pft_med_SQDU[4],range(-90,90,2), color=c2)
ax23.plot(pft_med_DQDU[4],range(-90,90,2), color=c3)
ax23.fill_betweenx(range(-90,90,2),pft_med_SQSU[4]-pft_std_SQSU[4], pft_med_SQSU[4]+pft_std_SQSU[4], facecolor='k', alpha=0.2)
ax23.fill_betweenx(range(-90,90,2),pft_med_SQDU[4]-pft_std_SQDU[4], pft_med_SQDU[4]+pft_std_SQDU[4], facecolor=c2, alpha=0.3)
ax23.fill_betweenx(range(-90,90,2),pft_med_DQDU[4]-pft_std_DQDU[4], pft_med_DQDU[4]+pft_std_DQDU[4], facecolor=c3, alpha=0.3)

line1, = ax24.plot(pft_med_SQSU[5],range(-90,90,2), color='k', label='Eppley')
line2, = ax24.plot(pft_med_SQDU[5],range(-90,90,2), color=c2, label='Kremer')
line3, = ax24.plot(pft_med_DQDU[5],range(-90,90,2), color=c3,label='Anderson')
leg = ax24.legend(handles = [line1,line2,line3],bbox_to_anchor=(.5, -1.2), loc="lower center", frameon=False)
for line in leg.get_lines():
    line.set_linewidth(3.0)
ax24 = plt.subplot(6,4,24)
ax24.set_xlabel('Zonal Mean\n(mg C $m^{-3}$)', fontsize=11)
ax24.fill_betweenx(range(-90,90,2),pft_med_SQSU[5]-pft_std_SQSU[5], pft_med_SQSU[5]+pft_std_SQSU[5], facecolor='k', alpha=0.2)
ax24.fill_betweenx(range(-90,90,2),pft_med_SQDU[5]-pft_std_SQDU[5], pft_med_SQDU[5]+pft_std_SQDU[5], facecolor=c2, alpha=0.3)
ax24.fill_betweenx(range(-90,90,2),pft_med_DQDU[5]-pft_std_DQDU[5], pft_med_DQDU[5]+pft_std_DQDU[5], facecolor=c3, alpha=0.3)

plt.gcf().text(0.34,-0.01, 'Mean Biomass over 0-240m', fontsize=12)
plt.gcf().text(0.17,0.91, '(a) Eppley', fontsize=14, weight='bold')
plt.gcf().text(0.4,0.91, '(b) Kremer', fontsize=14, weight='bold')
plt.gcf().text(0.61,0.91, '(c) Anderson', fontsize=14, weight='bold')
plt.gcf().text(0.84,0.91, '(d) ', fontsize=14, weight='bold')

#plt.savefig('figures/Figure4.pdf', bbox_inches='tight', transparent=True)


#%% Differences in biomass and PP (Figure S5)

biomass_SQSU_SQDU = np.divide(np.subtract(biomass_SQDU, biomass_SQSU),
                               biomass_SQSU)*100

biomass_SQSU_DQDU = np.divide(np.subtract(biomass_DQDU, biomass_SQSU),
                               biomass_SQSU)*100

biomass_SQDU_DQDU = np.divide(np.subtract(biomass_DQDU, biomass_SQDU),
                               biomass_SQDU)*100

cmap2 = matplotlib.cm.get_cmap('RdBu_r').copy()
cmap2.set_bad(color = 'grey', alpha = 1)

gridspec = dict(hspace=0,wspace=0.1,bottom=0.1, height_ratios=[1, 1, 1])
fig, [(ax1, ax4), (ax2, ax5), (ax3, ax6)] = plt.subplots(3, 2, figsize=(10, 9),sharey='row', sharex='col', gridspec_kw=gridspec)


ax1 = plt.subplot(3,2,1,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im1 = ax1.imshow(biomass_SQSU_SQDU, extent=(0,360,-90,90), cmap=cmap2, vmin=-50, vmax=50, origin='lower', transform=ccrs.PlateCarree())
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax1.set_yticks([-50,0,50])
ax1.set_yticklabels([50,0,-50])
ax1.set_title(r'$\bf{(a)}$'+'($C_{Kremer}$ - $C_{Eppley}$)/$C_{Eppley}$', loc='left')

ax2 = plt.subplot(3,2,3,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im2 = ax2.imshow(biomass_SQSU_DQDU, extent=(0,360,-90,90), cmap=cmap2, vmin=-50, vmax=50, origin='lower', transform=ccrs.PlateCarree())
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax2.set_yticks([-50,0,50])
ax2.set_yticklabels([50,0,-50])
ax2.set_title(r'$\bf{(c)}$'+'($C_{Anderson}$ - $C_{Eppley}$)/$C_{Eppley}$', loc='left')

ax3 = plt.subplot(3,2,5,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im3 = ax3.imshow(biomass_SQDU_DQDU, extent=(0,360,-90,90), cmap=cmap2, vmin=-50, vmax=50, origin='lower', transform=ccrs.PlateCarree())
fig.colorbar(im3,label='Phytoplankton Biomass Difference (%)', shrink=0.65, 
             orientation='horizontal', norm=colors.CenteredNorm(),cax=fig.add_axes([0.16,0.08,0.3,0.02]))
ax3.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax3.set_yticks([-50,0,50])
ax3.set_xlabel('Longitude')
ax3.set_yticklabels([50,0,-50])
ax3.set_title(r'$\bf{(e)}$'+'($C_{Anderson}$ - $C_{Kremer}$)/$C_{Kremer}$', loc='left')

ax4 = plt.subplot(3,2,2,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im4 = ax4.imshow(PP_SQSUc_SQDUc, extent=(0,360,-90,90), cmap=cmap2, vmin=-50, vmax=50, origin='lower', transform=ccrs.PlateCarree())
ax4.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax4.set_title(r'$\bf{(b)}$'+'($PP_{Kremer}$ - $PP_{Eppley}$)/$PP_{Eppley}$', loc='left')

ax5 = plt.subplot(3,2,4,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im5 = ax5.imshow(PP_SQSUc_DQDUc, extent=(0,360,-90,90), cmap=cmap2, vmin=-50, vmax=50, origin='lower', transform=ccrs.PlateCarree())
ax5.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
ax5.set_title(r'$\bf{(d)}$'+'($PP_{Anderson}$ - $PP_{Eppley}$)/$PP_{Eppley}$', loc='left')

ax6 = plt.subplot(3,2,6,projection=ccrs.PlateCarree(central_longitude=180, globe=None))
im6 = ax6.imshow(PP_SQDUc_DQDUc, extent=(0,360,-90,90), cmap=cmap2, vmin=-50, vmax=50, origin='lower', transform=ccrs.PlateCarree())
ax6.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', facecolor='k'))
fig.colorbar(im6,label='PP Variation (%)', shrink=0.65, 
             orientation='horizontal', norm=colors.CenteredNorm(),cax=fig.add_axes([0.56,0.08,0.3,0.02]))
ax6.set_xlabel('Longitude')
ax6.set_title(r'$\bf{(f)}$'+'($PP_{Anderson}$ - $PP_{Kremer}$)/$PP_{Kremer}$', loc='left')

fig.text(0.07,0.44, 'Latitude', va='center', rotation='vertical')

#plt.savefig('figures/FigureS5.pdf', bbox_inches='tight', transparent=True)

#%% Figure 2: Depth-Integrated biomass and primary production
# figure
conversion = 0.01201
bio_med_SQSU = np.nanmean(biomass_SQSU * conversion, axis=1) # units in mmol C / m2 * 1000 mmol/1 mol * 1 mol/12.01g * 1g/1000mg
bio_med_SQSU = bio_med_SQSU[11:87]
bio_std_SQSU = np.nanstd(biomass_SQSU * conversion, axis=1) 
bio_std_SQSU = bio_std_SQSU[11:87]
bio_med_SQDU = np.nanmean(biomass_SQDU * conversion, axis=1)
bio_med_SQDU = bio_med_SQDU[11:87]
bio_std_SQDU = np.nanstd(biomass_SQDU * conversion, axis=1)
bio_std_SQDU = bio_std_SQDU[11:87]
bio_med_DQDU = np.nanmean(biomass_DQDU * conversion, axis=1)
bio_med_DQDU = bio_med_DQDU[11:87]
bio_std_DQDU = np.nanstd(biomass_DQDU * conversion, axis=1)
bio_std_DQDU = bio_std_DQDU[11:87]

pp_med_SQSU = np.nanmean(PP_int_sum_SQSUc[8] * conversion * 3.154e+7, axis=1) # units in mmol C / m2 / s * 3.154e+7 seconds in a yr
pp_med_SQSU = pp_med_SQSU[11:87]
pp_std_SQSU = np.nanstd(PP_int_sum_SQSUc[8] * conversion * 3.154e+7, axis=1)
pp_std_SQSU = pp_std_SQSU[11:87]
pp_med_SQDU = np.nanmean(PP_int_sum_SQDUc[8] * conversion * 3.154e+7, axis=1)
pp_med_SQDU = pp_med_SQDU[11:87]
pp_std_SQDU = np.nanstd(PP_int_sum_SQDUc[8] * conversion * 3.154e+7, axis=1)
pp_std_SQDU = pp_std_SQDU[11:87] 
pp_med_DQDU = np.nanmean(PP_int_sum_DQDUc[8] * conversion * 3.154e+7, axis=1)
pp_med_DQDU = pp_med_DQDU[11:87]
pp_std_DQDU = np.nanstd(PP_int_sum_DQDUc[8] * conversion * 3.154e+7, axis=1)
pp_std_DQDU = pp_std_DQDU[11:87]

gridspec = dict(hspace=-0.5,wspace=0.1,bottom=0)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(7.5, 4),sharey='row', gridspec_kw=gridspec)
c2 = '#2A788EFF'
c3 = '#7AD151FF'

ax1.axhline(y=0, xmin=-0.5, xmax=3, color='grey', linestyle='--')
ax1.plot(bio_med_SQSU,range(-68,84,2), color='k')
ax1.fill_betweenx(range(-68,84,2),bio_med_SQSU-bio_std_SQSU, bio_med_SQSU+bio_std_SQSU, facecolor='k', alpha=0.1)
ax1.fill_betweenx(range(-68,84,2),bio_med_SQDU-bio_std_SQDU, bio_med_SQDU+bio_std_SQDU, facecolor=c2, alpha=0.1)
ax1.fill_betweenx(range(-68,84,2),bio_med_DQDU-bio_std_DQDU, bio_med_DQDU+bio_std_DQDU, facecolor=c3, alpha=0.1)
ax1.plot(bio_med_SQSU,range(-68,84,2), color='k')
ax1.plot(bio_med_SQDU,range(-68,84,2), color=c2)
ax1.plot(bio_med_DQDU,range(-68,84,2), color=c3)
ax1.set_xlabel('Biomass (mg C $m^{-2}$)', fontsize=11)
ax1.set_ylabel('Latitude', fontsize=11)

ax2.axhline(y=0, xmin=-50, xmax=200, color='grey', linestyle='--')
line1, = ax2.plot(pp_med_SQSU,range(-68,84,2), color='k',label='Eppley')
line2, = ax2.plot(pp_med_SQDU,range(-68,84,2), color=c2,label='Kremer')
line3, = ax2.plot(pp_med_DQDU,range(-68,84,2), color=c3,label='Anderson')
ax2.fill_betweenx(range(-68,84,2),pp_med_SQSU-pp_std_SQSU, pp_med_SQSU+pp_std_SQSU, facecolor='k', alpha=0.1)
ax2.fill_betweenx(range(-68,84,2),pp_med_SQDU-pp_std_SQDU, pp_med_SQDU+pp_std_SQDU, facecolor=c2, alpha=0.1)
ax2.fill_betweenx(range(-68,84,2),pp_med_DQDU-pp_std_DQDU, pp_med_DQDU+pp_std_DQDU, facecolor=c3, alpha=0.1)
leg = ax2.legend(handles = [line1,line2,line3],bbox_to_anchor=(.8, 0.77), loc="lower center", frameon=False)
for line in leg.get_lines():
    line.set_linewidth(3.0)
ax2.set_xlabel('Primary Production (mg C $m^{-2}$ $y^{-1}$)', fontsize=11)
fig.text(0.135,0.82, '(a)',fontweight='bold', fontsize=12)
fig.text(0.545,0.82, '(b)',fontweight='bold', fontsize=12)

#plt.savefig('figures/Figure2.pdf', bbox_inches='tight', transparent=True)
