import numpy as np
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"
})

# Constants
h = 6.62607015e-34  # Planck's constant (Joule second)
c = 3.0e8           # Speed of light in vacuum (m/s)
k_B = 1.380649e-23  # Boltzmann's constant (J/K)

# Planck function definition (in terms of frequency)
def planck_function_frequency(frequency, temperature):
    """
    Calculate the spectral radiance of a black body at a given frequency and temperature using the Planck function.
    
    Parameters:
    frequency: Frequency in Hz (numpy array or float)
    temperature: Temperature in Kelvin (float)
    
    Returns:
    Spectral radiance in W/m^2/sr/Hz
    """
    numerator = 2 * h * frequency**3 / c**2
    exponent = h * frequency / (k_B * temperature)
    denominator = np.exp(exponent) - 1
    return numerator / denominator
    
def tau_nu(nu,tau0,nu0,gamma):
    r = (nu-nu0)/gamma
    return tau0*np.exp(-r**2)

def observed_intensity(frequency,T1,T2,tau0,nu0,gamma):
    exponential = np.exp(-tau_nu(frequency,tau0,nu0,gamma))
    return planck_function_frequency(frequency,T1)*exponential+(1-exponential)*planck_function_frequency(frequency,T2)

# This is the intensity you'd get assuming the presence of a homogeneous medium

# Frequency range (in Hz)
frequencies = np.linspace(4e14, 6e14, 5000)  # from 1 THz to 1000 THz (infrared to visible range)

# Temperature in Kelvin
T1  = 6500
T2  = 5000
nu0 = 5e14
gamma = 2e13
tau0  = 10.

# Plot the result
f = plt.figure()
ax = f.add_subplot(111)
for tau0 in [1,3,10]:
    I     = observed_intensity(frequencies,T1,T2,tau0,nu0,gamma)
    ax.plot(frequencies * 1e-12, I, label=r'$I_{\nu,obs}$'+r'$\tau_0 = {}$'.format(tau0))
ax.plot(frequencies * 1e-12, planck_function_frequency(frequencies,T1), linestyle='dotted', label='background T = {} K'.format(T1))
ax.plot(frequencies * 1e-12, planck_function_frequency(frequencies,T2), linestyle='dotted', label='foreground T = {} K'.format(T2))
plt.legend()
plt.ylabel(f'Spectral Radiance (W m$-^2$sr$^-1$Hz$^-1$)')
plt.xlabel('Frequency (THz)')
ax = ax.twinx()
for tau0 in [1,3,10]:
    ax.plot(frequencies * 1e-12, tau_nu(frequencies,tau0,nu0,gamma), linestyle='dashed',
            label=r'$\tau_0 = {}$'.format(tau0))
plt.legend(loc='lower left')
plt.ylabel('optical depth')

plt.grid(True)
plt.savefig("img1.pdf", transparent = True)
    
# Temperature in Kelvin
T1  = 5000
T2  = 6500

# Plot the result
f = plt.figure()
ax = f.add_subplot(111)
for tau0 in [1,3,10]:
    I     = observed_intensity(frequencies,T1,T2,tau0,nu0,gamma)
    ax.plot(frequencies * 1e-12, I, label=r'$I_{\nu,obs}$'+r'$\tau_0 = {}$'.format(tau0))
ax.plot(frequencies * 1e-12, planck_function_frequency(frequencies,T1), linestyle='dotted', label='background T = {} K'.format(T1))
ax.plot(frequencies * 1e-12, planck_function_frequency(frequencies,T2), linestyle='dotted', label='foreground T = {} K'.format(T2))
plt.legend()
plt.ylabel('Spectral Radiance (W m$^-2$sr$^-1$ Hz$^-1$)')
plt.xlabel('Frequency (THz)')
ax = ax.twinx()
for tau0 in [1,3,10]:
    ax.plot(frequencies * 1e-12, tau_nu(frequencies,tau0,nu0,gamma), linestyle='dashed',
            label=r'$\tau_0 = {}$'.format(tau0))
plt.legend(loc='lower left')
plt.ylabel('optical depth')

plt.grid(True)
plt.savefig("img2.pdf", transparent = True)
plt.show()
