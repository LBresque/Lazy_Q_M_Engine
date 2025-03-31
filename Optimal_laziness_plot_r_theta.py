#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:02:18 2024

@author: Léa Bresque
"""
#Info 

#Code to plot the Optimal laziness probability of FIG. 4. in the main text and Fig. S8 in the SM

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
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
    #sol= (np.cos(theta)>=0)*sol_interm - (np.cos(theta)<0)*sol_interm
    #check dwdpl|pl=0 >0 
    partial_w_0 = -(1/ratio_loc)*(pj_loc*np.log(pj_loc)+(1-pj_loc)*np.log(1-pj_loc))/(np.sin(theta)*np.sin(theta/2)**2)
    partial_w_1 = -(4*np.sin(theta/2)**4/ratio_loc)*(pj_loc*np.log(pj_loc)+(1-pj_loc)*np.log(1-pj_loc))/(np.sin(theta)*np.sin(theta/2)**2)
    #return(1*(sol>=1)*(np.cos(theta)>0)+ 0*(sol<=0)*(np.cos(theta)>0) + 0*(sol>=1)*(np.cos(theta)<0)+ 1*(sol<=0)*(np.cos(theta)<0) + (partial_w_1<1)*(partial_w_0>1)*(sol>0)*(sol<1)*sol+ (partial_w_0>1)*(partial_w_1>1)*(sol>0)*(sol<1)*1+ (partial_w_0<1)*(partial_w_1<1)*(sol>0)*(sol<1)*0)
    return(1*(sol>=1)*(np.cos(theta)>0)+ 0*(sol<=0)*(np.cos(theta)>0) + (sol>0)*(sol<1)*sol)

def pl_opt_notrunc(theta,ratio_loc): # ratio = h w / kb T  (eq.(S97) in Supplementary Material)
    pj_loc = pj_fct(theta)
    inside= -ratio_loc*np.sin(theta)*np.sin(theta/2)**2/(pj_loc*np.log(pj_loc)+(1-pj_loc)*np.log(1-pj_loc))
    sol = (1-np.sqrt(inside))/np.cos(theta)
    return(sol)

    
'''

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
    
#%%%% cylindrical coordinate 
from matplotlib.colors import LogNorm

# Créez une grille de valeurs pour r et theta
r = np.linspace(0, 5, 100)  # Exemple de plage pour r
theta = np.linspace(0, 2*np.pi, 120)  # Plage complète pour theta

# Créez une grille 2D pour r et theta
R, Theta = np.meshgrid(r, theta)

# Calculez les valeurs de la fonction
Z = pl_opt(Theta,R) 
Z[Theta > np.pi] = 1

# Créez le plot en coordonnées polaires
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
c = ax.contourf(Theta, R, Z, cmap='viridis')
plt.colorbar(c)

# Affichez le plot
plt.show()

#%%%% 
'''

## theta in x axis 

#fig, ax = plt.subplots(figsize=(2.8,1)) #3.2,3.2*0.7
plt.gcf().set_dpi(600)

from matplotlib.colors import LinearSegmentedColormap

# Define the custom colormap
cmap_violet_white = LinearSegmentedColormap.from_list('violet_white', [(0, 'indigo'), (1, 'white')])
cmap_violet_white_red = LinearSegmentedColormap.from_list('violet_white_green', ["purple", "white", "#90EE90"])



# Créez une grille de valeurs pour r et theta
r = np.linspace(0, 3, 400)  # Exemple de plage pour r
theta = np.linspace(0, np.pi/2, 400)  # Plage complète pour theta

# Créez une grille 2D pour r et theta
#R, Theta = np.meshgrid(r, theta)
Theta, R = np.meshgrid(theta, r)

# Calculez les valeurs de la fonction
Z = pl_opt(Theta, R) 

vmin_loc= 0 
vmax_loc= 1
levels_loc = np.linspace(vmin_loc, vmax_loc, 11)
levels_loc_short = np.linspace(vmin_loc, vmax_loc, 11)
# Créez le plot 2D
plt.figure(figsize=(3.2,3.2*0.7))
contour = plt.contourf( Theta, R, Z, cmap=cmap_violet_white_red, levels= levels_loc, vmin= vmin_loc, vmax=vmax_loc, alpha=0.7) #bwr
plt.contour(Theta,R, Z,  levels=[l for l in levels_loc_short if l != 0], colors='black', alpha=1, linewidths=0.5) #'viridis'cmap_violet_white
cbar= plt.colorbar(contour)
cbar.set_ticks(np.linspace(vmin_loc, vmax_loc, 6)) #11
cbar.ax.tick_params(labelsize=9)
cbar.set_label('$p_l^{\\rm opt}$', fontsize=10)

line= plt.axhline(y=np.log(4),color='black', label=r'$\ln(4)$') #linewidth=1.5 linestyle='-'
#line.set_dashes([7, 3])
plt.plot(np.linspace(0,np.pi/2,10)[1:-1], np.ones(8)*np.log(4), color='black', linestyle=' ', label=r'$\ln(4)$', marker= 'x', markersize=4)

# Personnalisez les ticks pour l'axe theta
#theta_ticks = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi, 5 * np.pi / 4, 3 * np.pi / 2, 7 * np.pi / 4, 2 * np.pi]
#theta_labels = ['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$', r'$\frac{5\pi}{4}$', r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$', r'$2\pi$']
theta_ticks_short = [0, np.pi/8 , np.pi / 4, 3*np.pi/8, np.pi / 2]#[0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
theta_labels_short = ['0', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$']

#theta_labels_short = ['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$']

plt.xticks(theta_ticks_short, theta_labels_short)

# Personnalisez les ticks pour l'axe r
r_ticks = [0, 1, 2, 3] #np.log(4)
r_labels = ['0', '1', '2','3'] #r'$\ln(4)$'
plt.yticks(r_ticks, r_labels)

# Ajoutez des labels aux axes
plt.ylabel(r'$r$', fontsize=10)
plt.xlabel('$\\theta$', fontsize=10)
plt.yticks(fontsize=9)
plt.xticks(fontsize=9)
# Affichez le plot
#plt.show()

plt.savefig('pl_opt_r_theta_V18.pdf', format = 'pdf',bbox_inches='tight',pad_inches=0.025)



'''
#%%%%%%%%%
## theta in y axis 

fig, ax = plt.subplots(figsize=(3,3*0.7)) 
plt.gcf().set_dpi(600)

from matplotlib.colors import LinearSegmentedColormap

# Define the custom colormap
cmap_violet_white = LinearSegmentedColormap.from_list('violet_white', [(0, 'indigo'), (1, 'white')])
cmap_violet_white_red = LinearSegmentedColormap.from_list('violet_white_green', ["purple", "white", "#90EE90"])



# Créez une grille de valeurs pour r et theta
r = np.linspace(0, 3, 400)  # Exemple de plage pour r
theta = np.linspace(0, np.pi, 400)  # Plage complète pour theta

# Créez une grille 2D pour r et theta
R, Theta = np.meshgrid(r, theta)

# Calculez les valeurs de la fonction
Z = pl_opt(Theta, R) 

vmin_loc= 0 
vmax_loc= 1
levels_loc = np.linspace(vmin_loc, vmax_loc, 11)
# Créez le plot 2D
plt.figure(figsize=(8, 6))
contour = plt.contourf(R, Theta, Z, cmap=cmap_violet_white_red, levels= levels_loc, vmin= vmin_loc, vmax=vmax_loc, alpha=0.7) #bwr
plt.contour(R, Theta, Z,  levels=[l for l in levels_loc if l != 0], colors='black', alpha=1, linewidths=0.5) #'viridis'cmap_violet_white
cbar= plt.colorbar(contour)
cbar.set_ticks(np.linspace(vmin_loc, vmax_loc, 11))
cbar.ax.tick_params(labelsize=12)
cbar.set_label('$p_l^{\\rm opt}$', fontsize=12)

line= plt.axvline(x=np.log(4),color='black', linestyle='--', label=r'$\ln(4)$', linewidth=1.5)
line.set_dashes([7, 3])
plt.plot(np.ones(9)*np.log(4),np.linspace(0,np.pi,11)[1:-1], color='black', linestyle=' ', label=r'$\ln(4)$', marker= '+', markersize=8)

# Personnalisez les ticks pour l'axe theta
#theta_ticks = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi, 5 * np.pi / 4, 3 * np.pi / 2, 7 * np.pi / 4, 2 * np.pi]
#theta_labels = ['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$', r'$\frac{5\pi}{4}$', r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$', r'$2\pi$']
theta_ticks_short = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
theta_labels_short = ['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']

#theta_labels_short = ['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$']

plt.yticks(theta_ticks_short, theta_labels_short)

# Personnalisez les ticks pour l'axe r
r_ticks = [0, 1, np.log(4), 2, 3]
r_labels = ['0', '1',  r'$\ln(4)$','2','3']
plt.xticks(r_ticks, r_labels)

# Ajoutez des labels aux axes
plt.xlabel(r'$r$', fontsize=12)
plt.ylabel('$\\theta$', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
# Affichez le plot
#plt.show()

plt.savefig('pl_opt_r_theta_V6.pdf', format = 'pdf',bbox_inches='tight',pad_inches=0.025)
'''
