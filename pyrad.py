import pyradUtilities as utils
import pyradLineshape as ls
import pyradIntensity
import numpy as np
import matplotlib.pyplot as plt

c = 299792458.0
k = 1.38064852E-23
p0 = 1013.25
t0 = 296


def progressAlert():
    pass


def getCrossSection(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return obj.crossSection


def resetData(obj):
    obj.crossSection = np.copy(obj.yAxis)
    obj.progressCrossSection = False
    for child in obj:
        if not isinstance(child, Line):
            resetData(child)


def getAbsCoef(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return obj.absCoef


def getTransmissivity(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return obj.transmissivity


def getAbsorbance(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return obj.absorbance


def getEmissivity(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return obj.emissivity


def getGlobalIsotope(ID, isotopeDepth):
    globalIsoList = []
    for i in range(1, isotopeDepth + 1):
        globalIsoList.append(HITRAN_GLOBAL_ISO[ID][i])
    return globalIsoList


def printProgress(text, obj):
    layerName = obj.layer.name
    molName = obj.molecule.name
    isoName = obj.name
    print('Processing %s: %s; %s; isotope %s' % (text, layerName, molName, isoName))


def totalConcentration(layer):
    total = 0
    for molecule in layer:
        total += molecule.concentration
    return total


def totalLineList(obj):
    fullList = []
    if isinstance(obj, Isotope):
        return obj.linelist()
    for item in obj:
        fullList += totalLineList(item)
    return fullList


class Line:
    def __init__(self, wavenumber, intensity, einsteinA, airHalfWidth,
                 selfHalfWidth, lowerEnergy, tempExponent, pressureShift, parent):
        self.isotope = parent
        self.molecule = self.isotope.molecule
        self.layer = self.molecule.layer
        self.wavenumber = wavenumber
        self.intensity = intensity
        self.einsteinA = einsteinA
        self.airHalfWidth = airHalfWidth
        self.selfHalfWidth = selfHalfWidth
        self.lowerEnergy = lowerEnergy
        self.tempExponent = tempExponent
        self.pressureShift = pressureShift

    @property
    def broadenedLine(self):
        return self.wavenumber + self.pressureShift * self.layer.P / p0

    @property
    def lorentzHW(self):
        return ((1 - self.molecule.concentration) * self.airHalfWidth + self.molecule.concentration
                * self.selfHalfWidth) * (self.layer.P / p0) * (t0 / self.layer.T) ** self.tempExponent

    @property
    def gaussianHW(self):
        return self.broadenedLine * np.sqrt(2 * k * self.layer.T / self.isotope.molMass / c ** 2)


class Isotope(list):
    def __init__(self, number, molecule):
        super().__init__(self)
        params = utils.readMolParams(number)
        self.globalIsoNumber = params[0]
        self.shortName = params[1]
        self.name = 'Isotope %s' % self.globalIsoNumber
        self.molNum = params[2]
        self.isoN = params[3]
        self.abundance = params[4]
        self.q296 = params[5]
        self.gj = params[6]
        self.molMass = params[7]
        self.molecule = molecule
        self.layer = self.molecule.layer
        self.q = {}
        self.xAxis = np.copy(self.molecule.xAxis)
        self.yAxis = np.copy(self.molecule.yAxis)
        self.crossSection = np.copy(self.yAxis)
        self.progressCrossSection = False

    @property
    def P(self):
        return self.layer.P

    @property
    def T(self):
        return self.layer.T

    @property
    def depth(self):
        return self.layer.depth

    @property
    def rangeMin(self):
        return self.layer.rangeMin

    @property
    def rangeMax(self):
        return self.layer.rangeMax

    @property
    def resolution(self):
        return self.layer.resolution

    @property
    def distanceFromCenter(self):
        return self.layer.distanceFromCenter

    @property
    def absCoef(self):
        return self.crossSection * self.molecule.concentration * self.layer.P / 1E4 / k / self.layer.T

    @property
    def transmissivity(self):
        return np.exp(-self.absCoef * self.layer.depth)

    @property
    def emissivity(self):
        return 1 - self.transmissivity

    @property
    def absorbance(self):
        return np.log(1 / self.transmissivity)

    def getData(self):
        print('Getting data for %s, isotope #%s' % (self.molecule.name, self.globalIsoNumber))
        lineDict = utils.gatherData(self.globalIsoNumber, self.rangeMin, self.rangeMax)
        self.q = utils.getQData(self.globalIsoNumber)
        for line in lineDict:
            self.append(Line(line, lineDict[line]['intensity'], lineDict[line]['einsteinA'],
                             lineDict[line]['airHalfWidth'], lineDict[line]['selfHalfWidth'],
                             lineDict[line]['lowerEnergy'], lineDict[line]['tempExponent'],
                             lineDict[line]['pressureShift'], self))

    def createCrossSection(self):
        molecule = self.molecule
        layer = molecule.layer
        progress = 0
        i = 1
        alertInterval = int(len(self) / 20)
        crossSection = np.zeros(int((layer.rangeMax - layer.rangeMin) / layer.resolution))
        for line in self:
            print('Progress  <%s%s>\t nu=%s' % ('*' * i, '-' * (20 - i), line.wavenumber), end='\r', flush=True)
            if progress > i * alertInterval:
                i += 1
            progress += 1
            rightCurve = ls.pseudoVoigtShape(line.gaussianHW, line.lorentzHW,
                                             layer.resolution, layer.distanceFromCenter)
            intensity = pyradIntensity.intensityFactor(line.intensity, line.broadenedLine,
                                                       layer.T, line.lowerEnergy, self.q[layer.T], self.q296)
            arrayIndex = int((line.wavenumber - layer.rangeMin) / layer.resolution)
            arrayLength = len(crossSection) - 1
            if isBetween(arrayIndex, 0, arrayLength):
                crossSection[arrayIndex] = crossSection[arrayIndex] + rightCurve[0] * intensity
            for c in range(1, len(rightCurve) - 1):
                rightIndex = arrayIndex + c
                leftIndex = arrayIndex - c
                if isBetween(rightIndex, 0, arrayLength):
                    crossSection[rightIndex] += rightCurve[c] * intensity
                if isBetween(leftIndex, 0, arrayLength):
                    crossSection[leftIndex] += rightCurve[c] * intensity
        self.crossSection = crossSection
        self.progressCrossSection = True

    def linelist(self):
        lines = []
        for line in self:
            lines.append(line)
        return lines


class Molecule(list):
    def __init__(self, shortNameOrMolNum, layer, isotopeDepth=1, **abundance):
        super().__init__(self)
        self.layer = layer
        self.yAxis = np.copy(layer.yAxis)
        self.xAxis = layer.xAxis
        self.crossSection = np.copy(self.yAxis)
        try:
            int(shortNameOrMolNum)
            self.ID = int(shortNameOrMolNum)
            self.name = False
        except ValueError:
            self.name = shortNameOrMolNum
            self.ID = MOLECULE_ID[self.name]
        for isotope in getGlobalIsotope(self.ID, isotopeDepth):
            isoClass = Isotope(isotope, self)
            self.append(isoClass)
            if not self.name:
                self.name = isoClass.shortName
        self.concentration = 0
        self.progressCrossSection = False
        for key in abundance:
            if key == 'ppm':
                self.setPPM(abundance[key])
                self.concText = '%sppm' % abundance[key]
            elif key == 'ppb':
                self.setPPB(abundance[key])
                self.concText = '%sppb' % abundance[key]
            elif key == 'percentage':
                self.setPercentage(abundance[key])
                self.concText = '%s%%' % abundance[key]
            elif key == 'concentration':
                self.concentration = abundance[key]
                self.concText = '%s concentration' % abundance[key]
            else:
                print('Invalid concentration type. Use ppm, ppb, percentage, or concentration.')

    def __str__(self):
        return '%s: %s' % (self.name, self.concText)

    def setPercentage(self, percentage):
        self.concentration = percentage / 100

    def setPPM(self, ppm):
        self.concentration = ppm * 10**-6

    def setPPB(self, ppb):
        self.concentration = ppb * 10**-8

    def setConcentrationPercentage(self, percentage):
        self.setPPM(10000 * percentage)

    def getData(self):
        for isotope in self:
            isotope.getData()

    def createCrossSection(self):
        tempAxis = np.copy(self.yAxis)
        for isotope in self:
            tempAxis += getCrossSection(isotope)
        self.progressCrossSection = True
        self.crossSection = tempAxis

    @property
    def absCoef(self):
        return self.crossSection * self.concentration * self.layer.P / 1E4 / k / self.layer.T

    @property
    def transmissivity(self):
        return np.exp(-self.absCoef * self.layer.depth)

    @property
    def P(self):
        return self.layer.P

    @property
    def T(self):
        return self.layer.T

    @property
    def depth(self):
        return self.layer.depth

    @property
    def rangeMin(self):
        return self.layer.rangeMin

    @property
    def rangeMax(self):
        return self.layer.rangeMax

    @property
    def resolution(self):
        return self.layer.resolution

    @property
    def distanceFromCenter(self):
        return self.layer.distanceFromCenter


class Layer(list):
    layerList = []

    def __init__(self, depth, T, P, rangeMin, rangeMax, name=False, dynamicResolution=False):
        super().__init__(self)
        self.rangeMin = rangeMin
        self.rangeMax = rangeMax
        self.T = T
        self.P = P
        self.depth = depth
        self.distanceFromCenter = self.P / 1013.25 * 5
        if not dynamicResolution:
            self.resolution = utils.BASE_RESOLUTION
        else:
            self.resolution = 10**int(np.log10((self.P / 1013.25))) * .01
        Layer.layerList.append(self)
        self.xAxis = np.arange(rangeMin, rangeMax, self.resolution)
        self.yAxis = np.zeros(int((rangeMax - rangeMin) / self.resolution))
        self.crossSection = np.copy(self.yAxis)
        self.progressCrossSection = False
        if not name:
            name = 'layer %s' % Layer.layerList.index(self)
        self.name = name

    def __str__(self):
        return '%s; %s' % (self.name, '; '.join(str(m) for m in self))

    def createCrossSection(self):
        tempAxis = np.copy(self.yAxis)
        for molecule in self:
            tempAxis += getCrossSection(molecule)
        self.progressCrossSection = True
        self.crossSection = tempAxis

    @property
    def absCoef(self):
        tempAxis = np.copy(self.yAxis)
        for molecule in self:
            tempAxis += getAbsCoef(molecule)
        return tempAxis

    @property
    def transmissivity(self):
        return np.exp(-self.absCoef * self.depth)

    def resetData(self):
        for molecule in self:
            molecule.getData()

    def addMolecule(self, name, isotopeDepth=1, **abundance):
        molecule = Molecule(name, self, isotopeDepth, **abundance)
        self.append(molecule)
        if totalConcentration(self) > 1:
            print('**Warning : Concentrations exceed 1.')
        molecule.getData()
        return molecule


def returnPlot(obj, propertyToPlot):
    if propertyToPlot == "transmissivity":
        yAxis = getTransmissivity(obj), 1
    elif propertyToPlot == 'absorption coefficient':
        yAxis = getAbsCoef(obj), 0
    elif propertyToPlot == 'cross section':
        yAxis = getCrossSection(obj), 0
    elif propertyToPlot == 'absorbance':
        yAxis = getAbsorbance(obj), 0
    else:
        return False
    return yAxis


def isBetween(test, minValue, maxValue):
    if test >= minValue:
        if test <= maxValue:
            return True
    return False


def plot(obj, propertyToPlot, fill=True, individualColors=True):
    plt.figure(figsize=(10, 6), dpi=80)
    plt.subplot(111, facecolor='xkcd:dark grey')
    plt.xlabel('wavenumber cm-1')
    plt.margins(0.01)
    plt.subplots_adjust(left=.07, bottom=.08, right=.97, top=.90)
    plt.ylabel(propertyToPlot)
    plt.grid('grey', linewidth=.5, linestyle=':')
    plt.title('%s\nP: %smBars; T: %sK; depth: %scm' % (str(obj), obj.P, obj.T, obj.depth))
    yAxis, fillAxis = returnPlot(obj, propertyToPlot)
    fig, = plt.plot(obj.xAxis, yAxis, linewidth=1, color='w', alpha=.8, label=obj.name)
    plt.fill_between(obj.xAxis, fillAxis, yAxis, color='w', alpha=.3 * fill)
    handles = [fig]
    if type(yAxis) is bool:
        print('Invalid plot type. Choose "transmissivity", "absorption coefficient", "cross section", or "absorbance".')
        return False
    if isinstance(obj, list) and individualColors:
        if len(obj) > 6:
            print('More than 6 elements, only processing first 6...')
        for subPlot, color in zip(obj, COLOR_LIST):
            yAxis, fillAxis = returnPlot(subPlot, propertyToPlot)
            fig, = plt.plot(subPlot.xAxis, yAxis, linewidth=1, color=color, alpha=.8, label='%s' % subPlot.name)
            handles.append(fig)
            plt.fill_between(subPlot.xAxis, fillAxis, yAxis, color=color, alpha=.3 * fill)
    legend = plt.legend(handles=handles, frameon=False)
    text = legend.get_texts()
    plt.setp(text, color='w')
    plt.show()


HITRAN_GLOBAL_ISO = {1: {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 129},
                     2: {1: 7, 2: 8, 3: 9, 4: 10, 5: 11, 6: 12, 7: 13, 8: 14, 9: 121, 10: 15, 11: 120, 12: 122},
                     3: {1: 16, 2: 17, 3: 18, 4: 19, 5: 20},
                     4: {1: 21, 2: 22, 3: 23, 4: 24, 5: 25, },
                     5: {1: 26, 2: 27, 3: 28, 4: 29, 5: 30, 6: 31},
                     6: {1: 32, 2: 33, 3: 34, 4: 35},
                     7: {1: 36, 2: 37, 3: 38},
                     8: {1: 39, 2: 40, 3: 41},
                     9: {1: 42, 2: 43},
                     10: {1: 44},
                     11: {1: 45, 2: 46},
                     12: {1: 47, 2: 117},
                     13: {1: 48, 2: 49, 3: 50},
                     14: {1: 51, 2: 110},
                     15: {1: 52, 2: 53, 3: 107, 4: 108},
                     16: {1: 19, 2: 11, 3: 111, 4: 112},
                     17: {1: 56, 2: 113},
                     18: {1: 57, 2: 58},
                     19: {1: 59, 2: 60, 3: 61, 4: 62, 5: 63},
                     20: {1: 64, 2: 65, 3: 66},
                     21: {1: 67, 2: 68},
                     22: {1: 69, 2: 118},
                     23: {1: 70, 2: 71, 3: 72},
                     24: {1: 73, 2: 74},
                     25: {1: 75},
                     26: {1: 76, 2: 77, 3: 105},
                     27: {1: 78, 2: 106},
                     28: {1: 79},
                     29: {1: 80, 2: 119},
                     30: {1: 126},
                     31: {1: 81, 2: 82, 3: 83},
                     32: {1: 84},
                     33: {1: 85},
                     34: {1: 86},
                     35: {1: 127, 2: 128},
                     36: {1: 87},
                     37: {1: 88, 2: 89},
                     38: {1: 90, 2: 91},
                     39: {1: 92},
                     40: {1: 93, 2: 94},
                     41: {1: 95},
                     42: {1: 96},
                     43: {1: 116},
                     44: {1: 109},
                     45: {1: 103, 2: 115},
                     46: {1: 97, 2: 98, 3: 99, 4: 100},
                     47: {1: 114},
                     48: {1: 123},
                     49: {1: 124, 2: 125}}

COLOR_LIST = ['xkcd:bright orange',
              'xkcd:seafoam green',
              'xkcd:bright blue',
              'xkcd:salmon',
              'xkcd:light violet',
              'xkcd:green yellow']


MOLECULE_ID = {'h2o': 1, 'co2': 2, 'o3': 3, 'n2o': 4, 'co': 5,
               'ch4': 6, 'o2': 7, 'no': 8, 'so2': 9,
               'no2': 10, 'nh3': 11, 'hno3': 12, 'oh': 13,
               'hf': 14, 'hcl': 15, 'hbr': 16, 'hi': 17,
               'clo': 18, 'ocs': 19, 'h2co': 20, 'hocl': 21,
               'n2': 22, 'hcn': 23, 'ch3cl': 24, 'h2o2': 25,
               'c2h2': 26, 'c2h6': 27, 'ph3': 28, 'cof2': 29,
               'sf6': 30, 'h2s': 31, 'hcooh': 32, 'ho2': 33,
               'o': 34, 'clono2': 35, 'no+': 36, 'hobr': 37,
               'c2h4': 38, 'ch3oh': 39, 'ch3br': 40, 'ch3cn': 41,
               'cf4': 42, 'c4h2': 43, 'hc3n': 44, 'h2': 45,
               'cs': 46, 'so3': 47, 'c2n2': 48, 'cocl2': 49}
