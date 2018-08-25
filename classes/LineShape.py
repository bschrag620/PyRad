import scipy.constants as sc
import numpy as np
from matplotlib import pyplot as plt
from scipy.special import wofz
#import classes.PyradConstants as pr

t0 = 296
k = sc.k
c = sc.c
p0 = 1013.25

class LineShapeSetting():
    def __init__(self):
        self.wingWidth = 0

def gaussianHW(wavenumber, t, m):
    hW = wavenumber * np.sqrt(2 * k * t / m / c**2)
    return hW

def lorentzHW(airHalfWidth, selfHalfWidth, P, T, q, tExponent):

    hw = ((1-q) * airHalfWidth + q * selfHalfWidth) * (P / p0) * (t0 / T)**tExponent
    return hw

def gaussianLineShape(halfWidth, x):
    """Returns the right half of a gaussian curve, used for temp broadening in low pressure scenarios"""
    lineShape = np.sqrt(np.log(2) / np.pi) / halfWidth* np.exp(-(x / halfWidth)**2 * np.log(2))
    return lineShape

def lorentzLineShape(halfWidth, x):
    """Returns the right half of a lorentzian curve."""
    lineShape = halfWidth / (sc.pi * (x**2 + halfWidth**2))
    return lineShape


#   implemented
def pseudoVoigtShape(gHW, lHW, step):
    g = 2 *gHW
    l = 2 * lHW
    f = (g**5 + 2.69269 * g**4 * l + 2.42843 * g**3 * l**2 + 4.47163 * g**2 * l**3 + .07842 * g * l**4 + l**5)**(.2)
    n = 1.36603 * (l / f) - .47719 * (l / f)**2 + .11116 * (l / f)**3
    curveLength = np.arange(0, 4, step)
    gCurve = gaussianLineShape(f / 2, curveLength)
    lCurve = lorentzLineShape(f / 2, curveLength)
    psuedoVoigt = n * lCurve + (1 - n) * gCurve
    return psuedoVoigt


def broadenLineList(p, lineList):
    newList = {}
    for key in lineList:
        if key != 'Q':
            pressureShift = lineList[key]['pressureShift']
            newKey = key + pressureShift * p / p0
            newList[newKey] = lineList[key]    #transfers contents of old non-broadened dictionary to dictionary with a new pressure broadened key

    return newList, sorted(newList.keys())




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

"""voigt profile calculated from fadeeva function. used as a sanity check for the pseudo function as well as comparison for benchmarking.
borrowed from : https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/"""

def V(x, lHW, gHW):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    sigma = lHW / np.sqrt(2 * np.log(2))

    return np.real(wofz((x + 1j*gHW)/sigma/np.sqrt(2))) / sigma\
                                                           /np.sqrt(2*np.pi)

if __name__ == '__main__':
    """a bit of test code  for checking the voigt profile"""
    step = .001
    x = np.arange(-5, 5, step)
    gHW = .25
    lHW = .1
    fadeevaVoigt = V(x, gHW, lHW)
    pseudo = pseudoVoigtShape(gHW, lHW, step).tolist()
    pseudoAdjusted = np.zeros(len(fadeevaVoigt))
    position = int(len(pseudoAdjusted) / 2)
    while len(pseudo) > 0:
        pseudoAdjusted[position] = pseudo.pop(0)
        position += 1

    plt.plot(x, fadeevaVoigt, 'b', linewidth = .5)
    plt.plot(x, pseudoAdjusted, 'g', linewidth = .5)
    plt.show()