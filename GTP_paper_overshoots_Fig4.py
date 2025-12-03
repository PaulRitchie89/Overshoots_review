# -*- coding: utf-8 -*-
"""
Created on Tue May 27 13:13:23 2025

@author: pdlr201

Python script to create Figure 4 in "The implications of overshooting 1.5C
on Earth system tipping elements - a review"

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import rc
from matplotlib.patches import Patch
import matplotlib as mpl


from scipy.stats import exponnorm
from scipy.optimize import least_squares


    
def fun(x,TP_low,TP_cen,TP_high):
    return np.array([exponnorm.ppf(0.05,x[0],x[1],x[2])-TP_low, exponnorm.ppf(0.95,x[0],x[1],x[2])-TP_high,exponnorm.ppf(0.5,x[0],x[1],x[2])-TP_cen])
      
    

fontsize = 10
rc('font', **{'size' : fontsize, 'family' : 'serif'})
rc('text', usetex=True)

TP_threshs_cen = [1.5,1.5,1.8,3,3.5,4,7.5,1.2,2.8,4,1.5,2]
TP_threshs_low = [0.8,1,1.1,2,2,1.4,6,1,2,1.4,1,1.5]
TP_threshs_high = [3,3,3.8,6,6,8,10,1.5,3.5,5,2.3,3]
TP_tscales_cen = [10000,2000,10,2000,100,50,10000,10,50,100,200,200]
TP_tscales_low = [1000,500,5,500,50,15,10000,10,10,50,100,50]
TP_tscales_high = [15000,13000,50,10000,200,300,10000,10,500,100,300,1000]
TP_names = ['GrIS','WAIS','SPG','EASB','AF','AMOC','EAIS','CR','SAM','BFD','PFAT','GLCR']
TP_names_long = ['Greenland Ice Sheet','West Antarctic Ice Sheet','Subpolar Gyre','Marine Basins East Antarctica','Amazon Rainforest','AMOC','Non-marine East Antarctica','Warm-Water Coral Reefs','West African Monsoon','Boreal Forest','Land Permafrost','Mountain Glaciers']

idx=np.argsort(np.array(TP_tscales_cen))

TP_tscales_cen = np.array(TP_tscales_cen)[idx]
TP_threshs_cen = np.array(TP_threshs_cen)[idx]
TP_tscales_low = np.array(TP_tscales_low)[idx]
TP_threshs_low = np.array(TP_threshs_low)[idx]
TP_tscales_high = np.array(TP_tscales_high)[idx]
TP_threshs_high = np.array(TP_threshs_high)[idx]
TP_names = np.array(TP_names)[idx]
TP_names_long = np.array(TP_names_long)[idx]

sigma_threshs = np.zeros(len(TP_threshs_cen))
K_threshs = np.zeros(len(TP_threshs_cen))
mu_threshs = np.zeros(len(TP_threshs_cen))
sigma_tscales = np.zeros(len(TP_threshs_cen))
K_tscales = np.zeros(len(TP_threshs_cen))
mu_tscales = np.zeros(len(TP_threshs_cen))

np.random.seed(seed=1)

for l in range(len(TP_threshs_cen)):

    x0_threshs=[1,TP_threshs_cen[l],1]
    x0_tscales=[1,TP_tscales_cen[l],1]
    
    res_lsq_threshs = least_squares(fun, x0_threshs, args=(TP_threshs_low[l],TP_threshs_cen[l],TP_threshs_high[l]))
    res_lsq_tscales = least_squares(fun, x0_tscales, args=(TP_tscales_low[l],TP_tscales_cen[l],TP_tscales_high[l]))
    
    K_threshs[l], mu_threshs[l], sigma_threshs[l] = res_lsq_threshs.x[0], res_lsq_threshs.x[1], res_lsq_threshs.x[2]
    K_tscales[l], mu_tscales[l], sigma_tscales[l] = res_lsq_tscales.x[0], res_lsq_tscales.x[1], res_lsq_tscales.x[2]

# cmap = matplotlib.colormaps.get_cmap('RdYlBu')
cmap = plt.get_cmap('RdYlBu')

t_ovrs = [10,30,100,300,1000,3000,10000,30000,100000]
T_peaks = np.linspace(1,5,17)
dt = (np.log10(T_peaks[1])-np.log10(T_peaks[0]))#/2
dT = T_peaks[1]-T_peaks[0]
T_stab = 1.5

fig = plt.figure(figsize=(7.007874,10))#2.959))
# fig = plt.figure(figsize=(18,7.6))
ax = fig.add_subplot(311, projection = "3d")

colours = np.zeros((len(TP_names),4))

mpl.rcParams['hatch.linewidth'] = 0.1  # previous pdf hatch linewidth
hatches = ['','..','OO','xx']


for j in range(len(t_ovrs)):
    for k in range(len(T_peaks)):
        
        if T_peaks[k]>1.5:
        
            _zpos = 0
            
            for i in range(len(TP_names)):
                colours[i,:] = cmap(i/(len(TP_names)-1))
                
                if T_peaks[k] > TP_threshs_cen[i]:
                    
                    if TP_threshs_cen[i] < 1.5:
                        
                        if j == (len(t_ovrs)-1) and (k == len(T_peaks)-1):
                            bar = ax.bar3d(T_peaks[k]-dT/2, np.log10(t_ovrs[j])-dt/2, _zpos, dT, dt, 1, color=colours[i], edgecolor='k', linewidth=0.1, shade=False, label = TP_names[i], hatch=hatches[i%len(hatches)])
                            bar._facecolors2d  = bar._facecolor3d
                            bar._edgecolors2d  = bar._edgecolor3d
                        else:
                            ax.bar3d(T_peaks[k]-dT/2, np.log10(t_ovrs[j])-dt/2, _zpos, dT, dt, 1, color=colours[i], edgecolor='k', linewidth=0.1, shade=False, hatch=hatches[i%len(hatches)])
                        _zpos += 1    # add the height of each bar to know where to start the next
                    
                    else:
                        crit_t_ovr_squared = 6*(TP_tscales_cen[i]**2)*TP_threshs_cen[i]*(T_peaks[k] - T_stab)/((T_peaks[k] - TP_threshs_cen[i])**2)
                        
                        if t_ovrs[j] > np.sqrt(crit_t_ovr_squared):
                            if j == (len(t_ovrs)-1) and (k == len(T_peaks)-1):
                                bar = ax.bar3d(T_peaks[k]-dT/2, np.log10(t_ovrs[j])-dt/2, _zpos, dT, dt, 1, color=colours[i], edgecolor='k', linewidth=0.1, shade=False, label = TP_names[i], hatch=hatches[i%len(hatches)])
                                bar._facecolors2d  = bar._facecolor3d
                                bar._edgecolors2d  = bar._edgecolor3d
                            else:
                                ax.bar3d(T_peaks[k]-dT/2, np.log10(t_ovrs[j])-dt/2, _zpos, dT, dt, 1, color=colours[i], edgecolor='k', linewidth=0.1, shade=False, hatch=hatches[i%len(hatches)])
                            _zpos += 1    # add the height of each bar to know where to start the next



N = 1000
TP_thresh_rv = np.zeros((len(TP_names),N))
TP_tscale_rv = np.zeros((len(TP_names),N))

for i in range(len(TP_names)):        
    TP_thresh_rv[i,:] = exponnorm.rvs(K_threshs[i], mu_threshs[i], sigma_threshs[i], size=N)
    TP_tscale_rv[i,:] = exponnorm.rvs(K_tscales[i], mu_tscales[i], sigma_tscales[i], size=N)
    
    
    while (np.min(TP_thresh_rv[i,:]) < TP_threshs_low[i]) and (np.max(TP_thresh_rv[i,:]) > TP_threshs_high[i]):
        N2 = np.sum(TP_thresh_rv[i,:]<TP_threshs_low[i])
        TP_thresh_rv[i,:][TP_thresh_rv[i,:]<TP_threshs_low[i]] = exponnorm.rvs(K_threshs[i], mu_threshs[i], sigma_threshs[i], size=N2)
        N2 = np.sum(TP_thresh_rv[i,:]>TP_threshs_high[i])
        TP_thresh_rv[i,:][TP_thresh_rv[i,:]>TP_threshs_high[i]] = exponnorm.rvs(K_threshs[i], mu_threshs[i], sigma_threshs[i], size=N2)


    while (np.min(TP_tscale_rv[i,:]) < TP_tscales_low[i]) and (np.max(TP_tscale_rv[i,:]) > TP_tscales_high[i]):
        N2 = np.sum(TP_tscale_rv[i,:]<TP_tscales_low[i])
        TP_tscale_rv[i,:][TP_tscale_rv[i,:]<TP_tscales_low[i]] = exponnorm.rvs(K_tscales[i], mu_tscales[i], sigma_tscales[i], size=N2)
        N2 = np.sum(TP_tscale_rv[i,:]>TP_tscales_high[i])
        TP_tscale_rv[i,:][TP_tscale_rv[i,:]>TP_tscales_high[i]] = exponnorm.rvs(K_tscales[i], mu_tscales[i], sigma_tscales[i], size=N2)


ax2 = fig.add_subplot(312, projection = "3d")
  
for j in range(len(t_ovrs)):
    for k in range(len(T_peaks)):
        
        if T_peaks[k]>1.5:
        
            _zpos = 0
            
            for i in range(len(TP_names)):
                
                counter_ovr = 0
                
                for l in range(N):
                    
                    
                    if T_peaks[k] > TP_thresh_rv[i,l]:
                        
                        if TP_thresh_rv[i,l] < T_stab:
                            
                            counter_ovr += 1                            
                       
                        else:
                            crit_t_ovr_squared = 6*(TP_tscale_rv[i,l]**2)*TP_thresh_rv[i,l]*(T_peaks[k] - T_stab)/((T_peaks[k] - TP_thresh_rv[i,l])**2)
                            
                            if t_ovrs[j] > np.sqrt(crit_t_ovr_squared):
                                
                                counter_ovr += 1
                                
                    if counter_ovr > 100:
                        
                        if j == (len(t_ovrs)-1) and (k == len(T_peaks)-1):
                            bar = ax2.bar3d(T_peaks[k]-dT/2, np.log10(t_ovrs[j])-dt/2, _zpos, dT, dt, 1, color=colours[i], edgecolor='k', linewidth=0.1, shade=False, label = TP_names[i], hatch=hatches[i%len(hatches)])
                            bar._facecolors2d  = bar._facecolor3d
                            bar._edgecolors2d  = bar._edgecolor3d
                        else:
                            ax2.bar3d(T_peaks[k]-dT/2, np.log10(t_ovrs[j])-dt/2, _zpos, dT, dt, 1, color=colours[i], edgecolor='k', linewidth=0.1, shade=False, hatch=hatches[i%len(hatches)])
                        _zpos += 1    # add the height of each bar to know where to start the next
                        
                        break



ax3 = fig.add_subplot(313, projection = "3d")
  
for j in range(len(t_ovrs)):
    for k in range(len(T_peaks)):
        
        if T_peaks[k]>1.5:
        
            _zpos = 0
            
            for i in range(len(TP_names)):
                
                counter_ovr = 0
                
                for l in range(N):
                    
                    
                    if T_peaks[k] > TP_thresh_rv[i,l]:
                        
                        if TP_thresh_rv[i,l] < T_stab:
                            
                            counter_ovr += 1                            
                       
                        else:
                            crit_t_ovr_squared = 6*(TP_tscale_rv[i,l]**2)*TP_thresh_rv[i,l]*(T_peaks[k] - T_stab)/((T_peaks[k] - TP_thresh_rv[i,l])**2)
                            
                            if t_ovrs[j] > np.sqrt(crit_t_ovr_squared):
                                
                                counter_ovr += 1
                                
                    if counter_ovr > 10:
                        
                        if j == (len(t_ovrs)-1) and (k == len(T_peaks)-1):
                            bar = ax3.bar3d(T_peaks[k]-dT/2, np.log10(t_ovrs[j])-dt/2, _zpos, dT, dt, 1, color=colours[i], edgecolor='k', linewidth=0.1, shade=False, label = TP_names[i], hatch=hatches[i%len(hatches)])
                            bar._facecolors2d  = bar._facecolor3d
                            bar._edgecolors2d  = bar._edgecolor3d
                        else:
                            ax3.bar3d(T_peaks[k]-dT/2, np.log10(t_ovrs[j])-dt/2, _zpos, dT, dt, 1, color=colours[i], edgecolor='k', linewidth=0.1, shade=False, hatch=hatches[i%len(hatches)])
                        _zpos += 1    # add the height of each bar to know where to start the next
                        
                        break



def log_tick_formatter(val, pos=None):
    return f"10$^{{{int(val)}}}$" 

ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
ax2.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
ax2.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
ax3.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
ax3.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax3.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax3.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax3.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

ax.set_xlim(1.5,np.max(T_peaks))
ax.set_ylim(0.8,5.2)
ax.set_zlim(0,12)
ax2.set_xlim(1.5,np.max(T_peaks))
ax2.set_ylim(0.8,5.2)
ax2.set_zlim(0,12)
ax3.set_xlim(1.5,np.max(T_peaks))
ax3.set_ylim(0.8,5.2)
ax3.set_zlim(0,12)

ax.set_zticks([0,3,6,9,12])
ax2.set_zticks([0,3,6,9,12])
ax3.set_zticks([0,3,6,9,12])

ax.view_init(30,-125)
ax2.view_init(30,-125)
ax3.view_init(30,-125)

ax.tick_params(axis='x', labelsize=8, pad=-4)
ax.tick_params(axis='y', labelsize=8, pad=-4)
ax.tick_params(axis='z', labelsize=8, pad=-2)
ax2.tick_params(axis='x', labelsize=8, pad=-4)
ax2.tick_params(axis='y', labelsize=8, pad=-4)
ax2.tick_params(axis='z', labelsize=8, pad=-2)
ax3.tick_params(axis='x', labelsize=8, pad=-4)
ax3.tick_params(axis='y', labelsize=8, pad=-4)
ax3.tick_params(axis='z', labelsize=8, pad=-2)

ax.set_xlabel('Peak warming ($^oC$)',labelpad=-8,fontsize=fontsize)
ax2.set_xlabel('Peak warming ($^oC$)',labelpad=-8,fontsize=fontsize)
ax3.set_xlabel('Peak warming ($^oC$)',labelpad=-8,fontsize=fontsize)
ax.set_ylabel('Time over 1.5$^oC$\n(years)',labelpad=-4,fontsize=fontsize)
ax2.set_ylabel('Time over 1.5$^oC$\n(years)',labelpad=-4,fontsize=fontsize)
ax3.set_ylabel('Time over 1.5$^oC$\n(years)',labelpad=-4,fontsize=fontsize)
ax.zaxis.set_rotate_label(False)
ax2.zaxis.set_rotate_label(False)
ax3.zaxis.set_rotate_label(False)
ax.set_zlabel('Number of tipped\nelements',rotation=90,fontsize=fontsize,labelpad=-4)
ax2.set_zlabel('Number of tipped\nelements',rotation=90,fontsize=fontsize,labelpad=-4)
ax3.set_zlabel('Number of tipped\nelements',rotation=90,fontsize=fontsize,labelpad=-4)

fig.tight_layout()
# fig.subplots_adjust(bottom=0.3,top=1.0, left=0)
fig.subplots_adjust(bottom=0.0,top=1.0, left=0.0, right=0.8, wspace=0.0)

legend_elements = []
for j in range(len(TP_names_long)):
    legend_elements.append(Patch(facecolor=colours[j], edgecolor='k', linewidth=0.1, label=TP_names_long[j], hatch=hatches[j%len(hatches)]))


fig.legend(handles=legend_elements, ncol=1, loc='center', bbox_to_anchor=(0.83, 0.5), frameon=False, fontsize=fontsize, labelspacing=2)


fig.text(0.75,0.77,'\\textbf{Fast timescales}')
fig.text(0.75,0.22,'\\textbf{Slow timescales}')

fig.text(0.18,0.97,'\\textbf{Best estimate}')
fig.text(0.18,0.63,'\\textbf{10\% tipping risk}')
fig.text(0.18,0.29,'\\textbf{1\% tipping risk}')

fig.text(0.04,0.97,'\\textbf{(a)}')
fig.text(0.04,0.63,'\\textbf{(b)}')
fig.text(0.04,0.29,'\\textbf{(c)}')