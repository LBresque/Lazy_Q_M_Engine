#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:02:18 2024

@author: eriopse
"""
#Info 

#Code to plot the Optimal laziness probability of FIG. 4. in the main text and Fig. S8 in the SM

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=False)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

pl = 0.95
n0= 0 
nT = 10 
T=2 
kb = 1
n_cycles = 100
hw = 1
theta = np.pi/8 #the smaller : the smaller pj and the smaller the gain per cycle
Omega = 1 
cycles = np.arange(1,n_cycles+1, 1)
ratio = hw/(kb*T)

def pj_fct(theta):
    return(np.sin(theta/2)**2)
    
def pl_opt(theta,ratio_loc): # ratio = h w / kb T  (eq.(S97) in Supplementary Material)
    pj_loc = pj_fct(theta)
    inside= -ratio_loc*np.sin(theta)*np.sin(theta/2)**2/(pj_loc*np.log(pj_loc)+(1-pj_loc)*np.log(1-pj_loc))
    sol = (1-np.sqrt(inside))/np.cos(theta)
    return(1*(sol>=1)+ 0*(sol<=0) + (sol>0)*(sol<1)*sol)

    
theta_list = np.concatenate((np.arange(0,np.pi/2, np.pi/64)+np.pi/(128*10),np.array([np.pi/2-np.pi/(128*10)])))
ratios = np.array([0.5,1,2*np.log(2),1.5,3,10]) #ratio = h w / kb T 

#plot
plt.gcf().set_dpi(600)
marker_list= ['s','o', 'x', '+', 'v', '*','d']
line_list= ['dashed','dotted', 'solid', 'dashdot', (0,(1,1)), (0,(3.,1,1,1)),(0,(3,10,1,10))]
color_list = ['lightcoral', 'orangered', 'teal', 'olivedrab', 'cadetblue','sienna' ]

fig, ax = plt.subplots(figsize=(3.2,3.2*0.7)) 
m= 0
for r in ratios: 
   if r==2*np.log(2):
       ax.plot(theta_list,pl_opt(theta_list,np.log(4)), label='$r=\\ln(4)$', linestyle=line_list[m], marker=marker_list[m] ,color='black', clip_on=False, markevery=0.08,mfc='none', ms= 4) #marker=marker_list[m]
   else :
       ax.plot(theta_list,pl_opt(theta_list,r), linestyle=line_list[m],marker=marker_list[m], label='$r=%1.1f$' %r, clip_on=False, markevery=0.1, color=color_list[m], ms= 4)# mfc='none' 
   m+=1

plt.fill_between(theta_list, pl_opt(theta_list,np.log(4)), color= "powderblue",alpha= 0.2)
plt.fill_between(theta_list, pl_opt(theta_list,np.log(4)), y2= np.ones(33), color= "navajowhite",alpha= 0.2)
      
ax.set_xticks(np.arange(0,np.pi/2+ np.pi/8, np.pi/8))
ax.set_xticklabels(['$0$','$\\pi/8$', '$\\pi/4$', '$3\\pi/8$', '$\\pi/2$' ])   
ax.set_ylim([-.03,1.03])
ax.set_xlim([-0.03,np.pi/2+0.03])


ax.set_xlabel('$\\theta$', fontsize=12)
ax.set_ylabel('$p_l^{opt}$', fontsize=12)


ax.legend(
           loc='lower left',fontsize=9, ncol=3,
           markerfirst=False,handlelength=1.5,
           handletextpad=0.3,labelspacing=0.05,
           frameon=False, bbox_to_anchor=(-0.01, 0.98, 1.0, 0.102), #()
           columnspacing = 1)

    

#plt.savefig('fig_test_V1.pdf', format = 'pdf',bbox_inches='tight',pad_inches=0.025)



