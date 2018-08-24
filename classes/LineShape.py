import scipy.constants as sc
import numpy as np
import classes.PyradConstants as pr

t0 = pr.tRef
k = pr.k
c = pr.c
p0 = pr.pRef

class LineShapeSetting():
    def __init__(self):
        self.wingWidth = 0

def gaussianHW(wavenumber, t, m):
    hW = wavenumber * np.sqrt(2 * k * t / m / c**2)
    return hW

def lorentzHW(airHalfWidth, selfHalfWidth, P, T, q, tExponent):

    hw = ((1-q) * airHalfWidth + q * selfHalfWidth) * (P / p0) * (t0 / T)**tExponent
    return hw

def gaussianLineShape(halfWidth, step):
    """Returns the right half of a gaussian curve, used for temp broadening in low pressure scenarios"""
    x = 0
    y = np.sqrt(np.log(2) / np.pi) / halfWidth* np.exp(-(x / halfWidth)**2 * np.log(2))
    shape = []
    xRange = []
    tolerance = y / 500
    while y > tolerance:
        shape.append(y)
        xRange.append(x)
        x += step
        y = np.sqrt(np.log(2) / np.pi) / halfWidth* np.exp(-(x / halfWidth)**2 * np.log(2))
    return shape

def lorentzLineShape(halfWidth, step):
    """Returns the right half of a lorentzian curve."""
    x = 0
    y = halfWidth / (sc.pi * (x**2 + halfWidth**2))
    shape = []
    xRange = []

    #   change the tolerance to expand and contract how far the curve gets processed.
    #   Higher number = lower tolerance = more accuracy = more processing. Pick your poison.
    tolerance = y / 500

    while y > tolerance:
        shape.append(y)
        xRange.append(x)
        x += step
        y = (halfWidth / np.pi / (x ** 2 + halfWidth ** 2))
    return shape


#   to be implemented later.
def pseudoVoigtShape(values):

    return False

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