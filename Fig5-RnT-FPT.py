import numpy as np
#import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


#-----------------------------
# Define analytical functions
# ----------------------------
def SNR(pj,pl,nT,n0):
    num = (nT - n0) * (1.-pl) * (1.-pl+2.*pj*pl)**2
    den = 4.*pj*pl * (1.-pl+pj*pl) * (1.+pl-2.*pj*pl)
    return num /den

def SNRarr(pj,pl,nT,n0):
    arr = np.zeros(len(pj))
    for i in range(len(pj)):
        arr[i] = SNR(pj[i],pl,nT,n0)
    return arr

def fpt_rntplus_inz_feedback(z,n,n0,pj,pl):
    chi = (1.- pj) * pl / (1.- pj * pl)
    term = (1.+(1.-2.*pj)*pl*z**2.) / (2.*z*(1.-pj*pl)) 
    zp = term + np.sqrt( term**2.0 - chi ) 
    zm = term - np.sqrt( term**2.0 - chi ) 
    if (np.absolute(zp) < 1. and np.absolute(zm) > 1. ):
        zs = zp
    elif (np.absolute(zm) < 1. and np.absolute(zp) > 1. ):
        zs = zm    
    return np.power(zs,n-n0)*np.power(chi,n0-n) 

# Numerical inversion of Z transform
def myInvZ(tt, Err, func, *params):     
    rr = 10.0**(-Err/(2.0*tt))
    cf = 1.0/(2.0*tt*rr**tt)
    s = 0.+0.j
    for i in range(1,tt):
        z = rr*np.exp(1.0j*i*np.pi/tt)
        s += (-1.0)**i*np.real(func(*(z, *params)))
    s = 2.0*s + func(*(rr, *params)) + (-1.)**tt*func(*(-rr, *params))
    return np.real(cf*s)


# -----------------
# Parameter values
# -----------------
pj = 0.1    # probability of jump
xtgt = 4    # target site
x0 = 0      # initial condition
plarr = np.array([0.2,0.5,0.8,1.])
flarr = np.array(['RnT-FPT-nT4-n00-pj0-1-pl0-2-nt6-nh4e5.dat',
                  'RnT-FPT-nT4-n00-pj0-1-pl0-5-nt6-nh4e5.dat',
                  'RnT-FPT-nT4-n00-pj0-1-pl0-8-nt6-nh4e5.dat',
                  'RnT-FPT-nT4-n00-pj0-1-pl1-0-nt6-nh4e5.dat'])
clarr = np.array(['C0','C1','C2','C3'])
lsarr = np.array(['--','-.',':','-'])
mkarr = np.array(['o','v','s','^'])


# -----------------
# Plot
# -----------------
fig, ax1 = plt.subplots(figsize=(3,3*0.75))
ax2 = inset_axes(ax1, width="40%", height="40%")


# Plot SNR in the inset 
x = np.linspace(0.0001, 1, 5001)
ax2.plot(x, SNRarr(x,0.2,10,1),'C0-.', label='$p_l \!=\! 0.2$', lw = 1)
ax2.plot(x, SNRarr(x,0.5,10,1),'C1-',label='$ 0.5$', lw = 1)
ax2.plot(x, SNRarr(x,0.8,10,1),'C2--',label='$ 0.8$', lw = 1)

ax2.set_xlim([0,1])
ax2.set_ylim([0,40])
ax2.set_xticks([0,0.5,1],['$0$','$0.5$','$1$'],fontsize=7)
ax2.set_yticks([0,20,40])
ax2.set_xlabel('$p_j$',fontsize=7,labelpad=0.5)
ax2.set_ylabel('${\mathrm{SNR}}({\mathcal{T}})$',fontsize=7,labelpad=0.5)
ax2.tick_params(axis='both', which='major', labelsize=7)
ax2.legend(loc='upper right',handlelength=1.3,
           fontsize=7, markerfirst=False,
           handletextpad=0.5,labelspacing=0.1,
           frameon=False, bbox_to_anchor=(1.05, 1.08))



# Plot FPT in the main canvas
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([3,10**3])
ax1.set_ylim([1*10**(-5),1.])
ax1.set_xlabel(r'${\mathcal{T}}$',labelpad=0.1)
ax1.set_ylabel(r'$F_{{\mathcal{T}}}(n_T|n_0)$')
ax1.set_yticks([0.00001,0.0001,0.001,0.01,0.1,1])


ax1.bar(4,1,width=0.8,bottom=None,align='center',color='C4',alpha=0.6,label=r'$p_l\!=\!0.0$')

for idx in range(len(plarr)):
    pl = plarr[idx]
  
    f = open(flarr[idx], 'r')
    lines = f.readlines()
    x = np.array( [ np.float128(line.split()[0]) for line in lines if not line.startswith("#")] )
    y = np.array( [ np.float128(line.split()[1]) for line in lines if not line.startswith("#")] )
    f.close()
    
    ax1.plot(x,y,mkarr[idx],color=clarr[idx],ms=3,markeredgewidth=0.6,mfc='none',label=r'$p_l\!=\! %.1f$'%(plarr[idx]))
    
    xtt = np.arange(4,1000,2)
    count_tt = np.zeros(len(xtt))
    for i in range(len(xtt)):
        count_tt[i] = myInvZ(xtt[i], 10, fpt_rntplus_inz_feedback,xtgt,x0,pj,pl)
    ax1.plot(xtt,count_tt,'k',linestyle='--',linewidth=0.75)
    

ax1.legend(loc='upper right',fontsize=9, ncol=3,
           markerfirst=False,handlelength=1.2,
           handletextpad=0.3,labelspacing=0.1,
           frameon=False, bbox_to_anchor=(1.035, 1.22))

x=np.arange(300,900,5)
y=2.*np.power(x,-1.5)
ax1.plot(x,y,'-',color='Blue',lw=0.75)
ax1.text(500,0.0002,r'${\mathcal{T}}^{-\frac{3}{2}}$',fontsize=8)


plt.savefig('plot_RnT-FPT-with-inset.pdf', format='pdf', bbox_inches='tight',pad_inches=0.001)
