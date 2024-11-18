#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:07:55 2024

@author: Léa Bresque
"""
#Info
#code to obtain the red and blue contour plots in Fig S8 


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
          
plt.rc('text', usetex=False)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
 
#%%%%  Parameters
#----------- Change here the value of T and hw to obtain the different plots 
T=1 
kb = 1
hw = 1
ratio = hw/(kb*T)
#-----------

Omega = 1
n_cycles = 100
cycles = np.arange(1,n_cycles+1, 1)


marker_list= ['s','o', '+', 'x', 'v', '*','d']
line_list= ['solid','dotted', 'dashed', 'dashdot', (0,(1,1)), (0,(3.,1,1,1)),(0,(3,10,1,10))]
color_list = ['lightcoral', 'orangered', 'black', 'olivedrab', 'cadetblue','sienna' ]
m = 4#index in ratios = np.array([0.5,1,2*np.log(2),1.5,3,10])

#%%% Fonctions 

def Wtheta(theta): 
    return((hw/2)*np.sin(theta))

def pj_fct(theta):
    return(np.sin(theta/2)**2)

def n(t,pj,pl):
    n0= 0 
    fac_1 = (1-pl)*(t+1)/(1-pl+2*pl*pj)
    fac_2 = 2*pl*pj*(1-((1-2*pj)*pl)**(t+1))/((1-pl+2*pl*pj)**2)
    return(n0-1 +fac_1 + fac_2)


def dn(t, pj, pl):
    fac_1 = (1-pl)/(1-pl+2*pl*pj)
    #fac_2 = (2*pl*pj*(1-pl*(1-2*pj))*((1-2*pj)*pl)**t)/((1-pl+2*pl*pj)**2)
    fac_2 = 2*pl*pj*((1-2*pj)*pl)**t/(1-pl+2*pl*pj)
    return(fac_1 + fac_2)

def MSD(pl,pj,t): #<(n-n0)^2>  eq 75 of SM
    t1 = (t**2-1)*(1-pl)**4
    t2 = 4*pj*(t+1)*pl*((1-pl)**2*((1-pl)*(t-3)+3) + pj*pl*(1 -pl)*((1 - pl)*(t- 10) + 8) - 4*(pj*pl)**3 - 4*pj**2*(2-3*pl)*pl**2)
    t3 = -4*pj*pl*(2*((2*pj-1)*(1-pl)**2-(3*pj-2)*(1-pl)+pj)*(1-((1-2*pj)*pl)**(t+1)) - (1-pl)*(t+1)*(1-pl+2*pl*pj)*((1-2*pj)*pl)**(t+1))
    return(1+(1/(1-pl+2*pl*pj)**4)*(t1 + t2+t3))
    
def dnss(pj, pl):
    return((1-pl)/(1-pl+2*pl*pj))

def Wi(T, pj, pl): #per cycle
    return(-kb*T*(1-pl)*(pj*np.log(pj)+(1-pj)*np.log(1-pj)))
    #*np.log(2)

def Pss(theta,pl): #Omega,T,hw considered constants 
    pj_loc = pj_fct(theta)
    return((Omega/theta)*((hw/2)*np.sin(theta)*dnss(pj_loc,pl)-Wi(T,pj_loc,pl)))
    
def pl_opt(theta,ratio_loc): # ratio = h w / kb T
    pj_loc = pj_fct(theta)
    inside= -ratio_loc*np.sin(theta)*np.sin(theta/2)**2/(pj_loc*np.log(pj_loc)+(1-pj_loc)*np.log(1-pj_loc))
    sol = (1-np.sqrt(inside))/np.cos(theta)
    return(1*(sol>=1)+ 0*(sol<=0) + (sol>0)*(sol<1)*sol)

def darken_color(color, amount):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(amount *c[0], amount * c[1], amount *c[2])

    
pj = pj_fct(theta)
Wth = Wtheta(theta)

color_list = ['lightcoral', darken_color('orangered',0.8), 'black', 'olivedrab', 'cadetblue','sienna' ]

#%%% instantaneous average power 
n_pl = 141
pl_list = np.arange(0,n_pl, 1)/(n_pl-1) 
pl_list[-1]-= 0.001 #to avoid completely lazy
theta_list = np.arange(0,np.pi/2, np.pi/64)+np.pi/128
pj_list = pj_fct(theta_list)
Pss_table = np.zeros((n_pl,np.size(theta_list)))

i,j = 0,0
for th in theta_list: 
    i=0
    for pl in pl_list:
        Pss_table[i,j] = Pss(th,pl)
        i+= 1 #indice for pl
    j+=1

    
#fig, ax = plt.subplots()
#divnorm=colors.TwoSlopeNorm(vcenter=0.)
#cax = ax.imshow(Pss_table, cmap=plt.cm.seismic, norm=divnorm, extent=[theta_list[0], theta_list[-1],pl_list[0], pl_list[-1]], aspect=1, origin='lower') 
##ax.contour(Pss_table, norm=divnorm, extent=[theta_list[0], theta_list[-1],pl_list[0], pl_list[-1]], aspect=2, colors='white',alpha=0.5)
#cbar = fig.colorbar(cax, label='$\\langle P_{ss}\\rangle$') 
#ax.set_xlabel('$\\theta$')
#ax.set_ylabel('$p_l$')
#ax.plot(theta_list,pl_opt(theta_list,ratio), linestyle='-', color='black')
#plt.title('$\\Omega$ = %1.2f,  ' %Omega + '$\\hbar \omega$= %1.1f,  ' %hw + '$k_BT$ = %1.1f'  %T)
#plt.ylim([-0.02,1.02])


#%%% with coutour 

fig, ax = plt.subplots(figsize=(3.2,3.2*0.7)) 
divnorm=colors.TwoSlopeNorm(vcenter=0.)
cax = ax.contourf(Pss_table, cmap=plt.cm.seismic, norm=divnorm, extent=[theta_list[0], theta_list[-1],pl_list[0], pl_list[-1]])
cbar = fig.colorbar(cax, label='$\\langle P_{ss}\\rangle$') 
ax.set_xlabel('$\\theta$', fontsize= 12)
ax.set_ylabel('$p_l$', fontsize= 12)
ax.set_xticks(np.arange(0,np.pi/2+ np.pi/8, np.pi/8))
ax.set_xticklabels(['$0$','$\\pi/8$', '$\\pi/4$', '$3\\pi/8$', '$\\pi/2$' ])  
#plt.title('$\\Omega$ = %1.2f,  ' %Omega + '$\\hbar \omega$= %1.1f,  ' %hw + '$k_BT$ = %1.1f'  %T)
#plt.plot(theta_list[0:30],pl_list[index_max[0:30]], linewidth = 2, color='black')
#color_list[m]
plt.plot(theta_list,pl_opt(theta_list,ratio), color='black', linestyle=line_list[m],marker=marker_list[m], markevery=0.15)
plt.ylim([-0.01,1.02])
plt.xlim([-0.01,1.6])

#name = 'energetics_Omega%.2f' %Omega + '_hw%.2f' %hw + 'kT%.2f' %T + '.pdf'
#plt.savefig(name, dpi=1200, format = 'pdf',bbox_inches='tight',pad_inches=0.025)

