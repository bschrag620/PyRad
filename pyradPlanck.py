import numpy as np
np.seterr(divide='ignore', invalid='ignore')

k = 1.38064852E-23
c = 299792458.0
h = 6.62607004e-34
pi = 3.141592653589793
t0 = 296
p0 = 1013.25


def planck(a, b):
    # all plank functions follow the basic form a/(e**b - 1)

    return a / (np.exp(b) - 1)


def planckHz(Hz, temp):
    # planck blackbody intensity
    # input in in Hz, output is Wm-2sr-1Hz-1

    a = 2 * h * Hz**3 / c**2
    b = h * Hz / k / temp
    intensity = planck(a, b)
    return intensity


def planckWavelength(lam, temp):
    # planck blackbody intensity
    # input is wavelength in um, output is Wm-2sr-1um-1

    a = 2.0E24 * h * c ** 2 / (lam ** 5)
    b = 10 ** 6 * h * c / lam / k / temp
    intensity = planck(a, b)
    return intensity


def planckWavenumber(n, temp):
    # planck blackbody intensity
    # input is wavenumber in meters, output is wm-2sr-1(cm-1)-1
    a = 2E8 * h * c**2 * n**3
    b = 100 * h * c * n / k / float(temp)
    intensity = planck(a, b)
    return intensity


def reverseWavenumber(n, intensity):
    pass



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