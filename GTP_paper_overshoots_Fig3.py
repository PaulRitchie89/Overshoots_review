# -*- coding: utf-8 -*-
"""
Created on Tue May 27 13:48:36 2025

@author: pdlr201

Python script to create Figure 3 in "The implications of overshooting 1.5C
on Earth system tipping elements - a review"

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl
import matplotlib.colors as colors

from scipy.stats import exponnorm
from scipy.optimize import least_squares


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def gradientbars(bars,cmap):

      ax = bars[0].axes
      lim = ax.get_xlim()+ax.get_ylim()
      if cmap.name == 'hot':
          cmap = truncate_colormap(cmap, 0.2, 1.0)
      elif cmap.name == 'hot_r':
          cmap = truncate_colormap(cmap, 0.0, 0.8)
      for bar in bars:

          bar.set_zorder(1)
          bar.set_facecolor("none")
          x,y = bar.get_xy()
          w, h = bar.get_width(), bar.get_height()
          grad = np.atleast_2d(np.linspace(0,1,256))
          ax.imshow(grad.T, extent=[x,x+w,y,y+h], aspect="auto", zorder=0, norm=mpl.colors.NoNorm(vmin=0,vmax=1), cmap = cmap)
      ax.axis(lim)
     
def fun(x,TP_low,TP_cen,TP_high):
    return np.array([exponnorm.ppf(0.05,x[0],x[1],x[2])-TP_low, exponnorm.ppf(0.95,x[0],x[1],x[2])-TP_high,exponnorm.ppf(0.5,x[0],x[1],x[2])-TP_cen])
      
    

fontsize = 10
rc('font', **{'size' : fontsize, 'family' : 'serif'})
rc('text', usetex=True)

TP_threshs_cen = [1.5,1.5,1.8,3,3.5,4,7.5,1.2,2.8,4,2,1.5]
TP_threshs_low = [0.8,1,1.1,2,2,1.4,6,1,2,1.4,1.5,1]
TP_threshs_high = [3,3,3.8,6,6,8,10,1.5,3.5,5,3,2.3]
TP_tscales_cen = [10000,2000,10,2000,100,50,10000,10,50,100,200,200]
TP_tscales_low = [1000,500,5,500,50,15,10000,10,10,50,50,100]
TP_tscales_high = [15000,13000,50,10000,200,300,10000,10,500,100,1000,300]
TP_names = ['GrIS','WAIS','SPG','EASB','AF','AMOC','EAIS','CR','SAM','BFD','GLCR','PFAT']
TP_names_long = ['Greenland Ice Sheet','West Antarctic Ice Sheet','Subpolar Gyre','Marine Basins East Antarctica','Amazon Rainforest','AMOC','Non-marine East Antarctica','Warm-Water Coral Reefs','West African Monsoon','Boreal Forest','Mountain Glaciers','Land Permafrost']

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

for l in range(len(TP_threshs_cen)):

    x0_threshs=[1,TP_threshs_cen[l],1]
    x0_tscales=[1,TP_tscales_cen[l],1]
    
    res_lsq_threshs = least_squares(fun, x0_threshs, args=(TP_threshs_low[l],TP_threshs_cen[l],TP_threshs_high[l]))
    res_lsq_tscales = least_squares(fun, x0_tscales, args=(TP_tscales_low[l],TP_tscales_cen[l],TP_tscales_high[l]),ftol=1E-12)
    
    K_threshs[l], mu_threshs[l], sigma_threshs[l] = res_lsq_threshs.x[0], res_lsq_threshs.x[1], res_lsq_threshs.x[2]
    K_tscales[l], mu_tscales[l], sigma_tscales[l] = res_lsq_tscales.x[0], res_lsq_tscales.x[1], res_lsq_tscales.x[2]



t_ovrs = [10,30,100,300,1000,3000,10000,30000,100000]
T_peaks = np.linspace(1,5,17)
dt = (np.log10(T_peaks[1])-np.log10(T_peaks[0]))
dT = T_peaks[1]-T_peaks[0]
T_stab = 1.5

fig, ax = plt.subplots(1,2,sharex=True,figsize=(7.007874,4))



colours = np.zeros((len(TP_names),4))

ax[0].set_title('Global warming stabilises\nat peak warming', y=1.01, fontsize=fontsize,linespacing=1.5)
ax[1].set_title('Global warming exceeds 1.5$^oC$\nfor 100 years', y=1.01, fontsize=fontsize,linespacing=1.5)
ax[0].set_xlabel('Peak warming ($^oC$)')
ax[1].set_xlabel('Peak warming ($^oC$)')
ax[0].set_ylabel('Number of tipped elements')
ax[1].set_ylabel('Number of tipped elements')

ax[0].set_xlim(1, np.max(T_peaks))
ax[0].set_ylim(0,12)
ax[1].set_ylim(0,12)
ax[0].set_yticks(np.linspace(0,12,13))
ax[1].set_yticks(np.linspace(0,12,13))

t_ovrs2 = 100

N = 1000

counter_ovr = np.zeros((N,len(T_peaks)))
counter_stab = np.zeros((N,len(T_peaks)))
hist_ovr = np.zeros((len(TP_names)+1,len(T_peaks)))
hist_stab = np.zeros((len(TP_names)+1,len(T_peaks)))
TP_thresh_rv = np.zeros((len(TP_names),N))
TP_tscale_rv = np.zeros((len(TP_names),N))


levels = [0.000001,0.01,0.1,0.33,0.66,0.9,0.99,1.001]
cmap = plt.get_cmap('hot_r')
cmap = truncate_colormap(cmap, 0.1, 0.85)

# Take colors at regular intervals spanning the colormap.
colours = cmap(np.linspace(0, 1, len(levels)-1))

np.random.seed(seed=1)


for i in range(len(TP_names)):        
    TP_thresh_rv[i,:] = exponnorm.rvs(K_threshs[i], mu_threshs[i], sigma_threshs[i], size=N)
    TP_tscale_rv[i,:] = exponnorm.rvs(K_tscales[i], mu_tscales[i], sigma_tscales[i], size=N)
    
    count=0
    count2=0
    
    while (np.min(TP_thresh_rv[i,:]) < TP_threshs_low[i]) and (np.max(TP_thresh_rv[i,:]) > TP_threshs_high[i]):
        N2 = np.sum(TP_thresh_rv[i,:]<TP_threshs_low[i])
        TP_thresh_rv[i,:][TP_thresh_rv[i,:]<TP_threshs_low[i]] = exponnorm.rvs(K_threshs[i], mu_threshs[i], sigma_threshs[i], size=N2)
        N2 = np.sum(TP_thresh_rv[i,:]>TP_threshs_high[i])
        TP_thresh_rv[i,:][TP_thresh_rv[i,:]>TP_threshs_high[i]] = exponnorm.rvs(K_threshs[i], mu_threshs[i], sigma_threshs[i], size=N2)
        count += 1
        print(count)

    while (np.min(TP_tscale_rv[i,:]) < TP_tscales_low[i]) and (np.max(TP_tscale_rv[i,:]) > TP_tscales_high[i]):
        N2 = np.sum(TP_tscale_rv[i,:]<TP_tscales_low[i])
        TP_tscale_rv[i,:][TP_tscale_rv[i,:]<TP_tscales_low[i]] = exponnorm.rvs(K_tscales[i], mu_tscales[i], sigma_tscales[i], size=N2)
        N2 = np.sum(TP_tscale_rv[i,:]>TP_tscales_high[i])
        TP_tscale_rv[i,:][TP_tscale_rv[i,:]>TP_tscales_high[i]] = exponnorm.rvs(K_tscales[i], mu_tscales[i], sigma_tscales[i], size=N2)
        count2 += 1
        print(count2)
    

for k in range(len(T_peaks)):
    
    for j in range(N):
    
        for i in range(len(TP_names)):        
            
            if T_peaks[k] > TP_thresh_rv[i,j]:
                
                if TP_thresh_rv[i,j] <= T_stab:
                    
                    counter_ovr[j,k] += 1
                    
                else:
                
                    crit_t_ovr_squared_low = 6*(TP_tscale_rv[i,j]**2)*TP_thresh_rv[i,j]*(T_peaks[k] - T_stab)/((T_peaks[k] - TP_thresh_rv[i,j])**2)
                    
                    if t_ovrs2 > np.sqrt(crit_t_ovr_squared_low):
                        
                        counter_ovr[j,k] += 1
                
        counter_stab[j,k] = np.sum(TP_thresh_rv[:,j]<T_peaks[k])
  
    hist_stab[:,k],_ = np.histogram(counter_stab[:,k],bins=np.linspace(-0.5,len(TP_names)+0.5,len(TP_names)+2),density=True)
    hist_ovr[:,k],_ = np.histogram(counter_ovr[:,k],bins=np.linspace(-0.5,len(TP_names)+0.5,len(TP_names)+2),density=True)
    
    _zpos = 11.5
    count_stop = 0
    count_stop2 = 0
    
    for i in range(len(TP_names)+1):
        
        if (np.sum(hist_stab[len(TP_names)-i:,k]) > 0) and (count_stop==0):
            bar = ax[0].bar(T_peaks[k], 1, width=3*dT/4, bottom=_zpos, color=colours[np.sum(np.sum(hist_stab[len(TP_names)-i:,k])>levels)-1])#, edgecolor='k', linewidth=0.1, shade=False, label = TP_names[i], hatch=hatches[i%len(hatches)])
            if np.sum(hist_stab[len(TP_names)-i:,k]) > levels[-2]:
                count_stop = 1
        if (np.sum(hist_ovr[len(TP_names)-i:,k]) > 0) and (count_stop2==0) and (T_peaks[k] > 1.5):
            bar2 = ax[1].bar(T_peaks[k], 1, width=3*dT/4, bottom=_zpos, color=colours[np.sum(np.sum(hist_ovr[len(TP_names)-i:,k])>levels)-1])
            if np.sum(hist_ovr[len(TP_names)-i:,k]) > levels[-2]:
                count_stop2 = 1
        _zpos -= 1    # add the height of each bar to know where to start the next


import seaborn as sns
sns.despine()

fig.tight_layout()
fig.subplots_adjust(top=0.9,bottom=0.27,left=0.07,right=0.98,hspace=0.25,wspace=0.25)

fig.text(0.02,0.95,'\\textbf{(a)}')
fig.text(0.52,0.95,'\\textbf{(b)}')

ax4 = fig.add_axes([0.225, 0.1, 0.6, 0.03])
norm = mpl.colors.BoundaryNorm(levels, cmap.N)
cb2 = mpl.colorbar.ColorbarBase(ax4,cmap=cmap,norm=norm,extend='both',boundaries=levels,ticks=levels,orientation='horizontal')
cb2.set_label('Cumulative probability density')