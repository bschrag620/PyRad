import numpy as np

k = 1.38064852E-23
c = 299792458.0

h = 6.62607004e-34
pi = 3.141592653589793
t0 = 296
p0 = 1013.25

#   both the boltzmann equation and stimulated emissions rely on h*c/k .
#   Need to adjust it for units of cm, rather than m, hence the multiplication by 100
c2 = c * h * 100 / k


def boltzmannFactors(E, t):
    #   uses lower energy level, temperature, and c2 constant
    bF = np.exp(-c2 * E / t) / \
         np.exp(-c2 * E / t0)
    return bF


def stimulatedEmissions(wavenumber, t):
    #   uses wavenumber of the center line and t
    sE = (1 - np.exp(-c2 * wavenumber / t)) / \
         (1 - np.exp(-c2 * wavenumber / t0))
    return sE


def intensityFactor(intensity, wavenumber, t, lowerEnergy, q, q0):
    iF = intensity * (q0 / q) * (stimulatedEmissions(wavenumber, t)) * (boltzmannFactors(lowerEnergy, t))
    return iF
