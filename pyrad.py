import pyradUtilities as utils
import pyradLineshape as ls
import pyradIntensity
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

k = sc.k


def getGlobalIsotope(ID, isotopeDepth):
    globalIsoList = []
    for i in range(1, isotopeDepth + 1):
        globalIsoList.append(HITRAN_GLOBAL_ISO[ID][i])
    return globalIsoList


class Molecule:
    def __init__(self, name, ID, molecularWeight, isotopeDepth=1, **abundance):
        self.name = name
        self.ID = ID
        self.molecularWeight = molecularWeight
        self.isotopeDepth = isotopeDepth
        self.isoList = getGlobalIsotope(self.ID, self.isotopeDepth)
        self.info = {'linelist': {},
                     'q': {}}
        self.q = self.info['q']
        self.concentration = 0
        self.linelist = self.info['linelist']
        self.crossSection = np.array([])
        self.absCoefficient = np.array([])
        self.progressCrossSection = False
        for key in abundance:
            if key == 'ppm':
                self.setPPM(abundance[key])
            elif key == 'ppb':
                self.setPPB(abundance[key])
            elif key == 'percentage':
                self.setPercentage(abundance[key])
            elif key == 'concentration':
                self.concentration = abundance[key]
            else:
                print('Invalid concentration type. Use ppm, ppb, percentage, or concentration.')

    def setPercentage(self, perc):
        self.concentration = perc / 100

    def setPPM(self, ppm):
        #   use to set ppm, converts to percentage
        self.concentration = ppm * 10**-6

    def setPPB(self, ppb):
        #   use to set ppb, converts to percentage
        self.concentration = ppb * 10**-8

    def setConcentrationPercentage(self, perc):
        #   use to set concentration via a percentage, handy for wv
        self.setPPM(10000 * perc)

    def getData(self, layer):
        print('Getting data for %s' % self.name)
        for isotope in self.isoList:
            self.linelist.update(utils.gatherData(isotope, layer.rangeMin, layer.rangeMax))
            self.q[isotope] = utils.getQData(isotope)

    def createCrossSection(self, layer, distanceFromCenter=1):
        crossSection = np.zeros(int((layer.rangeMax - layer.rangeMin) / layer.resolution))
        orderedList = sorted(self.linelist.keys())
        for absoprtionLine in orderedList:
            lineDict = self.linelist[absoprtionLine]
            broadenedLine = ls.broadenLineList(layer.P, absoprtionLine, lineDict['pressureShift'])
            isotope = lineDict['isotope']
            globalIso = self.isoList[isotope - 1]
            qDict = self.q[globalIso]
            lhalfwidth = ls.lorentzHW(lineDict['airHalfWidth'], lineDict['selfHalfWidth'],
                                      layer.P, layer.T, self.concentration, lineDict['tempExponent'])
            ghalfwidth = ls.gaussianHW(broadenedLine, layer.T, self.molecularWeight)
            rightCurve = ls.pseudoVoigtShape(ghalfwidth, lhalfwidth, layer.resolution, distanceFromCenter)
            intensity = pyradIntensity.intensityFactor(lineDict['intensity'], broadenedLine, layer.T,
                                                       lineDict['lowerEnergy'],
                                                       qDict[layer.T], qDict[296])
            arrayIndex = int((absoprtionLine - layer.rangeMin) / layer.resolution)
            crossSection[arrayIndex] = crossSection[arrayIndex] + rightCurve[0] * intensity
            for c in range(1, len(rightCurve) - 1):
                rightIndex = arrayIndex + c
                leftIndex = arrayIndex - c
                if rightIndex < len(crossSection) - 1:
                    crossSection[rightIndex] += rightCurve[c] * intensity
                if leftIndex > 0:
                    crossSection[leftIndex] += rightCurve[c] * intensity
        self.crossSection = crossSection
        self.progressCrossSection = True
        return self.crossSection

    def createAbsorptionCoefficient(self, P, T):
        print('Creating absorption coefficient for %s' % self.name)
        self.absCoefficient = self.crossSection * self.concentration * P / 1E6 / k / T
        return self.absCoefficient

    def createTransmittance(self, x):
        print('Creating transmittance for %s' % self.name)
        transmittance = np.exp(-self.absCoefficient * x)
        return transmittance


def totalConcentration(molecules):
    total = 0
    for molecule in molecules:
        total += molecule.concentration
    return total


class Layer:
    layerList = []

    def __init__(self, depth, T, P, rangeMin, rangeMax, resolution=.01, name=False):
        self.rangeMin = rangeMin
        self.rangeMax = rangeMax
        self.resolution = resolution
        self.T = T
        self.P = P
        self.depth = depth
        self.layerComposition = []
        self.pressurePa = self.P * 100
        Layer.layerList.append(self)
        self.xAxis = np.arange(rangeMin, rangeMax, self.resolution)
        self.yAxis = np.zeros(int((rangeMax - rangeMin) / self.resolution))
        self.absCoef = self.yAxis
        self.transmittance = self.yAxis
        self.absorbance = self.yAxis
        self.progressAbsCoef = False
        self.progressTransmittance = False
        if not name:
            name = 'layer %s' % Layer.layerList.index(self)
        self.name = name

    def createCrossSection(self, distanceFromCenter=1):
        for molecule in self.layerComposition:
            print('Processing cross section for %s' % molecule.name)
            molecule.createCrossSection(self, distanceFromCenter)

    def getData(self):
        for molecule in self.layerComposition:
            molecule.getData(self)

    def addMolecules(self, *molecules):
        for molecule in molecules:
            self.layerComposition.append(molecule)
        if totalConcentration(molecules) > 1:
            print('**Warning : Concentrations exceed 1.')

    def createAbsorptionCoefficient(self, distanceFromCenter=1):
        for molecule in self.layerComposition:
            if not molecule.progressCrossSection:
                print('Absorption cross section for %s not yet processed, backtracking to crossSection...' % molecule.name)
                molecule.createCrossSection(self, distanceFromCenter)
        print('Creating absorption coefficient for %s' % self.name)
        for molecule in self.layerComposition:
            self.absCoef += molecule.createAbsorptionCoefficient(self.pressurePa, self.T)
        self.progressAbsCoef = True
        return self.absCoef

    def createTransmittance(self, distanceFromCenter=1):
        if not self.progressAbsCoef:
            print('Absorption coefficient not processed, backtracking to absorptionCoefficient')
            self.createAbsorptionCoefficient(distanceFromCenter)
        self.transmittance = np.exp(-self.absCoef * self.depth)
        self.progressTransmittance = True
        return self.transmittance

    def addMolecule(self, name, ID, isotopeDepth=1):
        molecule = Molecule(name, ID, isotopeDepth)
        self.layerComposition.append(molecule)
        if totalConcentration(self.layerComposition) > 1:
            print('**Warning : Concentrations exceed 1.')
        return molecule

    def createAbsorbance(self):
        if not self.progressTransmittance:
            print('Transmittance not processed, backtracking to absorptionCoefficient...')
            self.createTransmittance()
        self.absorbance = np.log(1 / self.transmittance)


def plot(xAxis, yAxis):
    plt.figure(figsize=(10,6), dpi=80)
    plt.subplot(111, facecolor='xkcd:dark grey')
    plt.xlabel('wavenumber cm-1')
    plt.margins(0.01)
    plt.subplots_adjust(left=.07, bottom=.08, right=.97, top=.97)
    plt.ylabel('transmittance')
    plt.grid('w', linewidth=.5)
    plt.plot(xAxis, yAxis, linewidth=1, color='orange')
    plt.show()

HITRAN_GLOBAL_ISO = {1: {1: 1,
                         2: 2,
                         3: 3,
                         4: 4,
                         5: 5,
                         6: 6,
                         7: 129},
                     2: {1: 7,
                         2: 8,
                         3: 9,
                         4: 10,
                         5: 11,
                         6: 12,
                         7: 13,
                         8: 14,
                         9: 121,
                         10: 15,
                         11: 120,
                         12: 122},
                     3: {1: 16,
                         2: 17,
                         3: 18,
                         4: 19,
                         5: 20},
                     4: {1: 21,
                         2: 22,
                         3: 23,
                         4: 24,
                         5: 25, },
                     5: {1: 26,
                         2: 27,
                         3: 28,
                         4: 29,
                         5: 30,
                         6: 31},
                     6: {1: 32,
                         2: 33,
                         3: 34,
                         4: 35},
                     7: {1: 36,
                         2: 37,
                         3: 38},
                     8: {1: 39,
                         2: 40,
                         3: 41},
                     9: {1: 42,
                         2: 43},
                     10: {1: 44},
                     11: {1: 45,
                          2: 46},
                     12: {1: 47,
                          2: 117},
                     13: {1: 48,
                          2: 49,
                          3: 50},
                     14: {1: 51,
                          2: 110},
                     15: {1: 52,
                          2: 53,
                          3: 107,
                          4: 108},
                     16: {1: 19,
                          2: 11,
                          3: 111,
                          4: 112},
                     17: {1: 56,
                          2: 113},
                     18: {1: 57,
                          2: 58},
                     19: {1: 59,
                          2: 60,
                          3: 61,
                          4: 62,
                          5: 63},
                     20: {1: 64,
                          2: 65,
                          3: 66},
                     21: {1: 67,
                          2: 68},
                     22: {1: 69,
                          2: 118},
                     23: {1: 70,
                          2: 71,
                          3: 72},
                     24: {1: 73,
                          2: 74},
                     25: {1: 75},
                     26: {1: 76,
                          2: 77,
                          3: 105},
                     27: {1: 78,
                          2: 106},
                     28: {1: 79},
                     29: {1: 80,
                          2: 119},
                     30: {1: 126},
                     31: {1: 81,
                          2: 82,
                          3: 83},
                     32: {1: 84},
                     33: {1: 85},
                     34: {1: 86},
                     35: {1: 127,
                          2: 128},
                     36: {1: 87},
                     37: {1: 88,
                          2: 89},
                     38: {1: 90,
                          2: 91},
                     39: {1: 92},
                     40: {1: 93,
                          2: 94},
                     41: {1: 95},
                     42: {1: 96},
                     43: {1: 116},
                     44: {1: 109},
                     45: {1: 103,
                          2: 115},
                     46: {1: 97,
                          2: 98,
                          3: 99,
                          4: 100},
                     47: {1: 114},
                     48: {1: 123},
                     49: {1: 124,
                          2: 125}}
