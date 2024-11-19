#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:22:34 2024

@author: Léa Bresque
"""

#%%%%% Information 

#Code to plot the Anomalous diffusion exponent

#Fig. 3 of the main text:  pj= 0.4
#Fig. S4 of the Supplementary Material:  pj= 0.2

#%%%%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
          
plt.rc('text', usetex=False)
#plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
 
#%%%%  Parameters

n_cycles = 10000 #number of cycles
cycles = np.arange(1,n_cycles+1, 1) 

#for plot lines
marker_list= ['s','o', '+', 'x', 'v', '*','d']
line_list= ['solid','dotted', 'dashed', 'dashdot', (0,(1,1)), (0,(3.,1,1,1)),(0,(3,10,1,10))]
color_list = ['lightcoral', 'orangered', 'black', 'olivedrab', 'cadetblue','sienna' ]

#%%% Mean square displacement function (eq 72 of SM)

def MSD(pl,pj,t): #<(nt-n0)^2>  (eq 72 of Supplementary Material)
    t1 = (t**2-1)*(1-pl)**4
    t2 = 4*pj*(t+1)*pl*((1-pl)**2*((1-pl)*(t-3)+3) + pj*pl*(1 -pl)*((1 - pl)*(t- 10) + 8) - 4*(pj*pl)**3 - 4*pj**2*(2-3*pl)*pl**2)
    t3 = -4*pj*pl*(2*((2*pj-1)*(1-pl)**2-(3*pj-2)*(1-pl)+pj)*(1-((1-2*pj)*pl)**(t+1)) - (1-pl)*(t+1)*(1-pl+2*pl*pj)*((1-2*pj)*pl)**(t+1))
    return(1+(1/(1-pl+2*pl*pj)**4)*(t1 + t2+t3))
    

#%%% Mean Standard deviation : diffusive, sub, super, ballistic 

fig, ax = plt.subplots(figsize=(3.2,3.2*0.5)) 
#------------------- Jump probability below to change depending on the plot
pj = 0.4 
#-------------------
m=-1 #maker index

for k in np.array([0,5,8,10]): 
    m+=1
    pl = k/10
    ax.plot(cycles, np.log(MSD(pl,pj,cycles))/np.log(cycles), label='$p_l=%1.2f$' %pl, linestyle=line_list[m])
    #Per_t2 = (1/(1-pl+2*pl*pj)**4)*((1-pl)**4+ 4*pj*pl*((1-pl)**3+pj*pl*(1-pl)**2))
    #plt.plot(cycles, np.log(MSD(pl,pj,cycles))/np.log(cycles)- np.log(Per_t2)/np.log(cycles), label='$p_l=%1.2f$' %pl)
    #plt.hlines(Per_t2,1,100)
    #plt.plot(cycles[0:20], MSD(pl,pj,cycles[0:20]), label='$p_l=%1.2f$' %pl)
#plt.plot(cycles,np.log(1+((1-2*pj)**(cycles+1)-1)/(2*pj**2)+(cycles+1)*(1-pj)/pj)/np.log(cycles), label='$p_l=1, th$', linestyle= ':')

cycles_short = np.array([2,10,100,1000,10000]) #just some cycle length values for the t and t^2 lines 
ax.plot(cycles_short,np.log(cycles_short)/np.log(cycles_short), label='$t$', linestyle= ':', marker= 'o', ms=5, color='gray')
ax.plot(cycles_short,np.log(cycles_short**2)/np.log(cycles_short), label='$t^2$', linestyle= ':', marker= 's', ms=5, color='black')   
ax.set_xscale('log')

ax.legend(loc='lower left',fontsize=9, ncol=3,
           markerfirst=False,handlelength=1.5,
           handletextpad=0.3,labelspacing=0.02,
           frameon=False, bbox_to_anchor=(-0.01, 0.95, 1.0, 0.102), #()
           columnspacing = 1)



ax.set_xlabel('Number of cycles $t$')
#ax.set_ylabel('$\\ln \\langle(n-n0)^2(t)\\rangle/\\ln t$')
ax.set_ylabel('$\\alpha(t)$')

    
#name = 'MSD_pj%.2f' %pj + '10000cycles.pdf'
#plt.savefig(name, dpi=1200, format = 'pdf',bbox_inches='tight',pad_inches=0.025)




    
