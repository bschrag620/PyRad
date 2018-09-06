import scipy.constants as sc
import numpy as np
from scipy.special import wofz
import pyradUtilities as utils


t0 = 296
k = sc.k
c = sc.c
p0 = 1013.25

cachedCurves = utils.getCurves('voigt')


def gaussianHW(wavenumber, t, m):
    hW = wavenumber * np.sqrt(2 * k * t / m / c**2)
    return hW


def lorentzHW(airHalfWidth, selfHalfWidth, P, T, q, tExponent):

    hw = ((1-q) * airHalfWidth + q * selfHalfWidth) * (P / p0) * (t0 / T)**tExponent
    return hw


def gaussianLineShape(halfWidth, xValue):
    """Returns the right half of a gaussian curve, used for temp broadening in low pressure scenarios"""
    lineShape = np.sqrt(np.log(2) / np.pi) / halfWidth * np.exp(-(xValue / halfWidth) ** 2 * np.log(2))
    return lineShape


def lorentzLineShape(halfWidth, xValue):
    """Returns the right half of a lorentzian curve."""
    lineShape = halfWidth / (sc.pi * (xValue**2 + halfWidth**2))
    return lineShape


def pseudoVoigtShape(gHW, lHW, dx, distanceFromCenter):
    gFW = 2 * gHW
    lFW = 2 * lHW
    fValue = (gFW**5 + 2.69269 * gFW**4 * lFW +
              2.42843 * gFW**3 * lFW**2 +
              4.47163 * gFW**2 * lFW**3 +
              .07842 * gFW * lFW**4 + lFW**5)**.2
    if fValue in cachedCurves:
        return cachedCurves[fValue]
    nValue = 1.36603 * (lFW / fValue) - .47719 * (lFW / fValue)**2 + .11116 * (lFW / fValue)**3
    curveLength = np.arange(0, distanceFromCenter, dx)
    gCurve = gaussianLineShape(fValue / 2, curveLength)
    lCurve = lorentzLineShape(fValue / 2, curveLength)
    psuedoVoigt = nValue * lCurve + (1 - nValue) * gCurve
    return psuedoVoigt


def broadenLineList(p, wavenumber, pressureShift):
    new = wavenumber + pressureShift * p / p0
    return new


def vvLineShape(halfwidth, waveCenter, step):
    x = waveCenter
    y = (halfwidth * waveCenter) / (sc.pi * waveCenter) * (1 / (waveCenter**2 + halfwidth**2) + 1/(waveCenter**2 + halfwidth**2))
    shape = []
    xRange = []
    tolerance = y / 500
    while y > tolerance:
        shape.append(y)
        xRange.append(x)
        x += step
        y = (halfwidth * x) / (sc.pi * waveCenter) * \
            (1 / ((x - waveCenter) ** 2 + halfwidth ** 2) + 1 / ((waveCenter + x) ** 2 + halfwidth ** 2))
    return shape


#   simply used to validate the shape of the pseudo curve in testing
def V(xValue, lHW, gHW):
    """Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha."""
    sigma = lHW / np.sqrt(2 * np.log(2))
    return np.real(wofz((xValue + 1j * gHW) / sigma / np.sqrt(2))) / sigma \
           / np.sqrt(2*np.pi)