from __future__ import print_function
import numpy as np
import pyradUtilities as utils


t0 = 296
k = 1.38064852E-23
c = 299792458.0
pi = 3.141592653589793
p0 = 1013.25

#cachedVoigt = utils.getCurves('voigt', utils.BASE_RESOLUTION)
cachedLorentz = utils.getCurves('lorentz', utils.BASE_RESOLUTION)
cachedGaussian = utils.getCurves('gaussian', utils.BASE_RESOLUTION)
print('\n', end='\r')


def gaussianHW(wavenumber, t, m):
    hW = wavenumber * np.sqrt(2 * k * t / m / c**2)
    return hW


def lorentzHW(airHalfWidth, selfHalfWidth, P, T, q, tExponent):
    hw = ((1-q) * airHalfWidth + q * selfHalfWidth) * (P / p0) * (t0 / T)**tExponent
    return hw


def gaussianLineShape(halfWidth, xValue):
    if isinstance(xValue, np.ndarray) and str(halfWidth) in cachedGaussian:
        length = len(xValue)
        cachedCurve = cachedGaussian[str(halfWidth)]
        if len(cachedCurve) >= length:
            return cachedCurve[:length]
    """Returns the right half of a gaussian curve, used for temp broadening in low pressure scenarios"""
    lineShape = np.sqrt(np.log(2) / np.pi) / halfWidth * np.exp(-(xValue / halfWidth) ** 2 * np.log(2))
    cachedGaussian[halfWidth] = lineShape
    return lineShape


def lorentzLineShape(halfWidth, xValue):
    if isinstance(xValue, np.ndarray) and str(halfWidth) in cachedLorentz:
        length = len(xValue)
        cachedCurve = cachedLorentz[str(halfWidth)]
        if len(cachedCurve) >= length:
            return cachedCurve[:length]
    """Returns the right half of a lorentzian curve."""
    lineShape = halfWidth / (pi * (xValue**2 + halfWidth**2))
    cachedLorentz[halfWidth] = lineShape
    return lineShape


def pseudoVoigtShape(gHW, lHW, dx, distanceFromCenter):
    gFW = 2 * gHW
    lFW = 2 * lHW
#    cacheKey = '%s:%s' % (gFW, lFW)
#    length = int(distanceFromCenter / dx)
#    if cacheKey in cachedVoigt:
#        cachedCurve = cachedVoigt[cacheKey]
#        if len(cachedCurve) >= length:
#            return cachedCurve[:length]
    fValue = (gFW**5 + 2.69269 * gFW**4 * lFW +
              2.42843 * gFW**3 * lFW**2 +
              4.47163 * gFW**2 * lFW**3 +
              .07842 * gFW * lFW**4 + lFW**5)**.2
    nValue = 1.36603 * (lFW / fValue) - .47719 * (lFW / fValue)**2 + .11116 * (lFW / fValue)**3
    curveLength = np.arange(0, distanceFromCenter, dx)
    gCurve = gaussianLineShape(fValue / 2, curveLength)
    lCurve = lorentzLineShape(fValue / 2, curveLength)
    pseudoVoigt = nValue * lCurve + (1 - nValue) * gCurve
#    cachedVoigt[cacheKey] = pseudoVoigt
    return pseudoVoigt


def broadenLineList(p, wavenumber, pressureShift):
    new = wavenumber + pressureShift * p / p0
    return new


def vvLineShape(halfwidth, waveCenter, step):
    x = waveCenter
    y = (halfwidth * waveCenter) / (pi * waveCenter) * (1 / (waveCenter**2 + halfwidth**2) + 1/(waveCenter**2 + halfwidth**2))
    shape = []
    xRange = []
    tolerance = y / 500
    while y > tolerance:
        shape.append(y)
        xRange.append(x)
        x += step
        y = (halfwidth * x) / (pi * waveCenter) * \
            (1 / ((x - waveCenter) ** 2 + halfwidth ** 2) + 1 / ((waveCenter + x) ** 2 + halfwidth ** 2))
    return shape


def writeCacheToFile():
    print('Writing line shapes to file.')
    utils.writeCurveToFile(cachedGaussian, 'gaussian', utils.BASE_RESOLUTION)
    utils.writeCurveToFile(cachedLorentz, 'lorentz', utils.BASE_RESOLUTION)
#   utils.writeCurveToFile(cachedVoigt, 'voigt', utils.BASE_RESOLUTION)

#   simply used to validate the shape of the pseudo curve in testing

'''def V(xValue, lHW, gHW):
    """Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha."""
    sigma = lHW / np.sqrt(2 * np.log(2))
    return np.real(wofz((xValue + 1j * gHW) / sigma / np.sqrt(2))) / sigma \
           / np.sqrt(2*np.pi)'''
