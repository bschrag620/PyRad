import scipy.constants as spc
import matplotlib.pyplot as plt
import numpy as np

##initilize the plot with some basic parameters:
#plt.grid(True)
#plt.margins(x=0, y=.01)

## some basic constants
c = spc.c   #speed of light in m
h = spc.h   #planck's constant
k = spc.k   #boltzmann constant

        #all plank functions follow the basic form a/(e**b - 1)
def planck(a,b):
    return a / (np.exp(b) - 1)

    #planck blackbody intensity
    #input in in Hz, output is Wm-2sr-1Hz-1
def planckHz(Hz, T):
    a = 2 * h * Hz**3 / c**2
    b = h * Hz / k / T
    intensity = planck(a,b)
    return intensity

    #planck blackbody intensity
    #input is wavelength in um, output is Wm-2sr-1um-1
def planckWavelength(lam, T):
    a = 2.0E24 * h * c ** 2 / (lam ** 5)
    b = 10**6 * h * c / lam / k / T
    intensity = planck(a,b)
    return intensity

    # planck blackbody intensity
    # input is wavenumber in meters, output is wm-2sr-1(cm-1)-1
def planckWavenumber(n, T):
    a = 2E8 * h * c**2 * n**3
    b = 100 * h * c * n / k / T
    intensity = planck(a, b)
    return intensity

def spectrumIntensity(self):
    intensity = planckWavenumber(self.vRange, self.t)
    return intensity

def spectrumWavelength(self):
    return self.vRange

def roundTo(self, n):
    self.vRange *= 10**n
    self.vRange = self.vRange.astype(int)
    self.vRange = self.vRange.astype(float)
    self.vRange /= 10**n

def length(self):
    return len(self.spectrumWavelength())




# Example usage and plot
##create arange array of wavenumbers 1-250000, units meters
#wavenumbers = np.arange(1, 2500, 1)

##calculate intensities for wavenumbers
#intensity = planckWavenumber(wavenumbers, T)

##add intensities to the plt
#plt.plot(wavenumbers, intensity, 'r', linewidth = .5)

##transmitted intensities calculated
#transmitted = intensity * .8

#plt.plot(wavenumbers, transmitted, 'b', linewidth = .5)


#plt.show()