import numpy as np 
import matplotlib.pyplot as plt 
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"
})
from matplotlib import cm
from astropy import constants as const

m_p = const.m_p.value #kg
c = const.c.value
h = const.h.value
m_e = const.m_e.value #kg
k = const.k_B.value #J/K
T0 = 1.5e6 #K
a = 1e9 #m
N0 = 3e7 #cm^-3
G = const.G.value #dimensioni giuste
M_s = const.M_sun.value
#l = G*M_s*m/(2*k*T0*a)
#a = a*1e2 #conversione a cm

def max(v, T):
    return (m_e/(2 * np.pi * k * T))**(3/2) * np.exp(-m_e * v**2 / (2 * k * T))*4*np.pi*v**2

def B(f, T):
    return ((2*h*f**3)/c**2)*(np.exp(h*f/(k*T))-1)**(-1)

def Bl(l, T):
    return B(c/l, T)*((c/l)**2)/c

v = np.linspace(0, 5e7, 1000) #m/s
#l = c/f
T = np.array([10, 1e2, 1e3, 1e4, 1e5, 1e6]) #array temperature
colors = plt.cm.viridis(np.linspace(0, 1, len(T)))  # color map

'''for idx,i in enumerate(T):
    plt.grid(which="both", ls="dashed", color="grey")
    plt.plot(v/c*1e2, f(v,i), label = f'Distribuzione di Maxwell-Boltzmann per $T$ = {i/1e6:.1f} $\cdot 10^6$ K', color = colors[idx])
    plt.xlabel(f'Fattore Beta di Lorentz [\%]')
    plt.ylabel('Distribuzione di probabilità')
    plt.legend(loc = 'upper right')
#plt.show()
plt.savefig("electronmax.pdf", transparent = True)'''

for idx,i in enumerate(T):
    plt.grid(which="both", ls="dashed", color="grey")
    if i == 10 or i == 100:
        f = np.linspace(1e8, 1e14, 1000)
    elif i == 1e4 or i==1e3:
        f = np.linspace(1e8, 1e16, 1000)
    else:
        f = np.linspace(1e8, 1e18, 1000)
    plt.plot(f, B(f, i), label = f'$T$ = {i:.1e} K', color = colors[idx])
    plt.xlabel(f'Frequency [Hz]')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-25, 1e12)
    plt.ylabel(f'Specific intensity [erg s$^2$ m$^-2$ Hz$^-2$]')
    plt.legend(loc = 'upper left')
#plt.show()
plt.savefig("blackbody.pdf", transparent = True)