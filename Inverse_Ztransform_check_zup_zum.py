import numpy as np
import scipy
import matplotlib.pyplot as plt
import time
from matplotlib import colors

# pi is the probability to use information processing, i.e., pi= 1- pl 

npi = 21
list_pi= np.linspace(0,10,npi)/10

#function definitions 

#First passage time generating function 
def FPT(u,nT,n0,pj,pi): 
    chi = (1.-pj)*(1.-pi)/(1.-pj+pj*pi)
    tt = (1.+(1.-2.*pj)*(1.-pi)*u**2.)/(2.*u*(1.-pj+pj*pi))
    zm = tt - np.sqrt( tt**2.0 - chi )
    zp = tt + np.sqrt( tt**2.0 - chi )
    z_sol = zm*(np.abs(zp)>np.abs(zm)) + zp*(np.abs(zp)<np.abs(zm)) 
    return(chi**(n0-nT)*z_sol**(nT-n0))

#roots z_u^+ and z_u^- given in Eq (41) or the SM
def zp_zm(u,nT,n0,pj,pi):
    chi = (1.-pj)*(1.-pi)/(1.-pj+pj*pi)
    tt = (1.+(1.-2.*pj)*(1.-pi)*u**2.)/(2.*u*(1.-pj+pj*pi))
    zm = tt - np.sqrt( tt**2.0 - chi )
    zp = tt + np.sqrt( tt**2.0 - chi )
    return(zp, zm)

#numerical method to obtain the time distribution from the generating function (see ref 32 of main text)
def myInvu(t, Err, func, *params):
    rr = 10.0**(-Err/(2.0*t))
    cf = 1.0/(2.0*t*rr**t)
    s = 0.+0.j
    for i in range(1,t):
        #print('it', i)
        u = rr*np.exp(1.0j*i*np.pi/t)
        s += (-1.0)**i*np.real(func(*(u, *params)))
    s = 2.0*s + func(*(rr, *params)) + (-1.)**t*func(*(-rr, *params))
    return np.real(cf*s)

#%% plot FPT theory with respect to time 
nbr_times=600
n0=1
nT=7
pj=0.1 
pi= 0.2

tt = np.arange(1,nbr_times+1,1)
FPT_list = np.zeros(nbr_times)

for idx in range(nbr_times-1):
    FPT_list[idx] = myInvu(tt[idx], 3, FPT,nT,n0,pj,pi)

#plt.plot(np.log(tt),np.log(FPT_list), linestyle= '',  marker = '+')
plt.plot(tt,FPT_list, linestyle= '',  marker = '+')
plt.ylabel('Probability FPT, log scale')
plt.xlabel('Time, log scale')
#plt.ylim([-5,0.1])
plt.ylim([0,1.2])
plt.xlim([0,15])

#%%% plot FPT theory with respect to time  for different values of pi, i.e., of pl 
marker_list = ['o', '*', 'D', 'x', 'h', 'H', 's', '8', '+','v', '<' ]

i = 0
for pi_loc in list_pi: 
    for idx in range(nbr_times-1):
        FPT_list[idx] = myInvu(tt[idx], 3, FPT,nT,n0,pj,pi_loc)
    plt.plot(tt,FPT_list, linestyle= '', label= '$p_l= %1.1f$' %(1-pi_loc), marker = marker_list[i])
    i+=1 


plt.title('$p_j =$ %.1f' %pj)
plt.ylabel('FPT Probability')
plt.xlabel('Time (nbr of cycles)')
plt.xscale('log')
plt.yscale('log')
plt.ylim([10**(-7),1])
plt.xlim([0,1000])
plt.legend()


#%%% color map of the modulus of z_u^+ and z_u^-, as well as of chi, for different modulus of u and a random phase, different pi, i.e., pl and a fixed pj

nu = 100
list_u =(np.linspace(-10,10,nu)/10)*np.exp(1j*np.random.random(1)*2*np.pi) #*np.exp(1j*np.random.random(nu)*2*np.pi) #np.linspace(-10,10,nu)/10 # #np.linspace(-10,10,npi)/10
zp_list = np.zeros((npi,nu)) + 1j*np.zeros((npi,nu))
zm_list =  np.zeros((npi,nu)) + 1j*np.zeros((npi,nu))

pj=0.1
u_var = 0 

for u_loc in list_u: 
    pi_var=0 
    for pi_loc in list_pi: 
        zp_list[pi_var][u_var] , zm_list[pi_var][u_var] = zp_zm(u_loc,nT,n0,pj,pi_loc)
        pi_var+=1
    u_var+=1

plt.subplot(1, 3, 1)
plt.imshow(np.abs(zp_list),vmin=0, vmax=2, cmap='seismic', aspect='auto', extent=[-1,1,1,0])
plt.title('$|z_p|$')
plt.xlabel('u')
plt.ylabel('$p_i$')
plt.gca().get_yaxis().set_visible(True)

plt.subplot(1, 3, 2)
plt.imshow(np.abs(zm_list),vmin=0, vmax=2, cmap='seismic', aspect='auto', extent=[-1,1,1,0])
plt.title('$|z_m|$')
plt.xlabel('u')
#cbar = plt.colorbar()
#cbar.set_label('Modulus', rotation=90)
#plt.suptitle('$p_j=$%1.2f ; ' %pj + '$n_T=$%1.2f ; ' %nT +'$n_0=$%1.2f' %n0)
#plt.show()
plt.gca().get_yaxis().set_visible(False)

plt.subplot(1, 3, 3) #np.abs(zp_list)*np.abs(zm_list)  #((np.abs(zp_list)<=1) & (np.abs(zm_list)>=1)) | ((np.abs(zp_list)>=1) & (np.abs(zm_list)<=1))
plt.imshow(np.abs(zp_list)*np.abs(zm_list),vmin=0, vmax=2, cmap='seismic', aspect='auto', extent=[-1,1,1,0])
plt.title('$|z_m|.|z_p|$')
plt.xlabel('u')
cbar = plt.colorbar()
cbar.set_label('Modulus', rotation=90)
plt.suptitle('$p_j=$%1.2f ; ' %pj + '$n_T=$%1.2f ; ' %nT +'$n_0=$%1.2f' %n0)
plt.gca().get_yaxis().set_visible(False)
plt.show()

#%%% modulus of z_u^+ and z_u^- with respect to the real part of u
#list_u = np.random.random(nu)*np.exp(1j*np.random.random(nu)*2*np.pi)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)

plt.subplot(1, 2, 1)
for i in range(1):
    plt.plot(np.real(list_u), np.abs(zp_list[i][:]), linestyle='none', marker= '+')
plt.xlabel('Real(u)')
plt.ylabel('$|z_p|$')
plt.hlines(1, -1,1, colors='black')
plt.vlines(0, 0,5, colors='black')
plt.ylim([-0.3,5])

plt.subplot(1, 2, 2)
for i in range(1):
    plt.plot(np.real(list_u), np.abs(zm_list[i][:]), linestyle='none', marker= '+')
plt.xlabel('Real(u)')
plt.ylabel('$|z_m|$')
plt.hlines(1, -1,1, colors='black')
plt.vlines(0, 0,5, colors='black')
plt.ylim([-0.3,5])

#%%% modulus of z_u^+ and z_u^- with respect to the imaginary part of u
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)

plt.subplot(1, 2, 1)
for i in range(11):
    plt.plot(np.imag(list_u), np.abs(zp_list[i][:]), linestyle='none', marker= '+')
plt.xlabel('Im(u)')
plt.ylabel('$|z_p|$')
plt.hlines(1, -1,1, colors='black')
plt.vlines(0, 0,5, colors='black')
plt.ylim([-0.3,10])

plt.subplot(1, 2, 2)
for i in range(11):
    plt.plot(np.imag(list_u), np.abs(zm_list[i][:]), linestyle='none', marker= '+')
plt.xlabel('Im(u)')
plt.ylabel('$|z_m|$')
plt.hlines(1, -1,1, colors='black')
plt.vlines(0, 0,5, colors='black')
plt.ylim([-0.3,10])

#%%  How tt evolves 
def tt(): 
    return()

def f_z(tt, chi): 
    return(tt * np.sqrt(tt**2 - chi))

list_fz = np.zeros((201,201))+ 1j*np.zeros((201,201))
list_tt = np.linspace(-5,5,201)
list_chi = np.linspace(0,10,201)

i_chi = 0 
for chi in list_chi: 
    i_tt = 0
    for tt in list_tt: 
        list_fz[i_tt][i_chi] = f_z(tt,chi)
        i_tt+=1
    i_chi+=1
    
plt.imshow(np.abs(list_fz), vmin=0, vmax = 2, cmap= colors.ListedColormap(['black', 'red']), extent=[0,10,-10,10])
plt.plot(list_chi[20],list_tt[100], 'r+')
cbar = plt.colorbar()
cbar.set_label('Modulus', rotation=90)
