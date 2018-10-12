import os
import pyradUtilities as utils
import pyradLineshape as ls
import pyradIntensity
import pyradPlanck
import numpy as np
import matplotlib.pyplot as plt
# import pyradInteractive


c = 299792458.0
k = 1.38064852E-23
h = 6.62607004e-34
pi = 3.141592653589793
R = 8.3144598
t0 = 296
p0 = 1013.25
avo = 6.022140857E23


def wvPressure(temperature, RH):
    wvP = RH * wvSaturationVaporPressure(temperature) / 100
    return wvP


def wvRH2AH(temperature, pressure, RH):
    wvP = wvPressure(temperature, RH)

    absHum = wvP / 461.5 * temperature
    print(wvP, absHum)
    return absHum


def wvMixingRatio(temperature, pressure, RH):
    e = wvPressure(temperature, RH)
    r = .622 * e / (pressure - e)
    return r


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def wvSaturationVaporPressure(temperature):
    temperature -= 273
    tC = 647.096
    pC = 220640
    a1 = -7.85951783
    a2 = 1.84408259
    a3 = -11.7866497
    a4 = 22.6807411
    a5 = -15.9618719
    a6 = 1.80122502
    # v = 1 - temperature / tC
    # satVP = pC * np.exp(tC / temperature * (a1 * v + a2 * v**1.5 + a3 * v**3 + a4 * v**3.5 + a5 * v**4 + a6 * v**7.5))
    satVP = 6.11 * np.exp(17.3 * temperature / (temperature + 237.3))
    return satVP


def wvConcentrationFromRH(layer, RH):
    temperature = layer.T
    pressure = layer.P
    concentration = wvMixingRatio(temperature, pressure, RH)
    return concentration

def linear(baseValue, baseHeight, rate, height):
    newValue = baseValue + rate * (height - baseHeight)
    return newValue


def integrateSpectrum(spectrum, unitAngle=pi, res=utils.BASE_RESOLUTION):
    value = np.sum(np.nan_to_num(spectrum))
    value = value * unitAngle * res
    return value


def getCrossSection(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return obj.crossSection


def resetCrossSection(obj):
    if not obj.progressCrossSection:
        return
    obj.crossSection = np.zeros(int((obj.rangeMax - obj.rangeMin) / utils.BASE_RESOLUTION))
    obj.progressCrossSection = False
    for child in obj:
        if not isinstance(child, Line):
            resetCrossSection(child)


def resetData(obj):
    # clears the existing line data from parent object down to isotope, and then reloads the data using getData
    # use this if layer ranges get changed. Will also clear the cross section data of the obj.
    for child in obj:
        if isinstance(child, Isotope):
            while len(child) > 0:
                child.pop()
            child.getData()
        else:
            resetData(child)
    resetCrossSection(obj)


def getAbsCoef(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return obj.absCoef


def getTransmittance(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return obj.transmittance


def getOpticalDepth(obj):
    if not obj.progressCrossSection:
        obj.createCrossSection()
    return -np.log(obj.transmittance)


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


def convertLength(value, units):
    if units == 'cm':
        return value
    elif units in ['m', 'meter']:
        return value * 100
    elif units in ['ft', 'feet']:
        return value * 30.48
    elif units in ['in', 'inch']:
        return value * 2.54


def convertPressure(value, units):
    if units == 'mbar':
        return value
    elif units in ['atm', 'atmospheres', 'atmosphere']:
        return value * 1013.25
    elif units in ['b', 'bar']:
        return value * 1000
    elif units in ['pa', 'pascal', 'pascals']:
        return value / 100


def convertRange(value, units):
    if units == 'cm-1':
        return value
    elif units in ['um', 'micrometers', 'micrometer']:
        return 10000 / value


def convertTemperature(value, units):
    if units[0].upper() == 'K':
        return value
    elif units[0].upper() == 'C':
        return value + 273
    elif units[0].upper() == 'F':
        return (value - 32) * 5 / 9 + 273


def interpolateArray(hiResXAxis, loResXAxis, loResYValues):
    hiResY = np.interp(hiResXAxis, loResXAxis, loResYValues)
    return hiResY


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
        return (float((1 - self.molecule.concentration) * self.airHalfWidth + self.molecule.concentration
                * self.selfHalfWidth) * (self.layer.P / p0) * (t0 / self.layer.T) ** self.tempExponent)

    @property
    def gaussianHW(self):
        return self.broadenedLine * np.sqrt(2 * k * self.layer.T / self.isotope.molMass / c ** 2)


class Isotope(list):
    def __init__(self, number, molecule):
        super(Isotope, self).__init__(self)
        params = utils.readMolParams(number)
        self.globalIsoNumber = params[0]
        self.shortName = params[1]
        self.name = 'Isotope %s' % self.globalIsoNumber
        self.molNum = params[2]
        self.isoN = params[3]
        self.abundance = params[4]
        self.q296 = params[5]
        self.gj = params[6]
        self.molmass = params[7]
        self.molecule = molecule
        self.layer = self.molecule.layer
        self.q = {}
        self.crossSection = np.copy(self.layer.crossSection)
        self.lineSurvey = np.zeros(int((self.layer.rangeMax - self.layer.rangeMin) / utils.BASE_RESOLUTION))
        self.progressCrossSection = False

    @property
    def P(self):
        return self.layer.P

    @property
    def molMass(self):
        return self.molmass / 1000 / avo

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
    def transmittance(self):
        return np.exp(-self.absCoef * self.layer.depth)

    @property
    def emissivity(self):
        return 1 - self.transmittance

    @property
    def emittance(self):
        return self.emissivity

    @property
    def absorbance(self):
        return np.log10(1 / self.transmittance)

    @property
    def yAxis(self):
        return np.copy(self.layer.yAxis)

    @property
    def xAxis(self):
        return np.copy(self.layer.xAxis)

    def getData(self, lineSurvey=False, verbose=True):
        if verbose:
            print('Getting data for %s, %s, isotope %s' % (self.layer.name, self.molecule.name, self.globalIsoNumber), end='\r', flush=True)
        lineDict = utils.gatherData(self.globalIsoNumber, self.layer.effectiveRangeMin, self.layer.effectiveRangeMax)
        self.q = utils.getQData(self.globalIsoNumber)
        for line in lineDict:
            # need to soft code the intensity value, depending on resolution setting
            if lineDict[line]['intensity'] > 1E-22:
                self.append(Line(line, lineDict[line]['intensity'], lineDict[line]['einsteinA'],
                             lineDict[line]['airHalfWidth'], lineDict[line]['selfHalfWidth'],
                             lineDict[line]['lowerEnergy'], lineDict[line]['tempExponent'],
                             lineDict[line]['pressureShift'], self))
        # self.createLineSurvey()

    def createCrossSection(self):
        molecule = self.molecule
        layer = molecule.layer
        progress = 0
        i = 1
        alertInterval = int(len(self) / 20)
        crossSection = np.copy(self.yAxis)
        trackGauss = 0
        trackLorentz = 0
        trackVoigt = 0
        if len(self) == 0:
            self.progressCrossSection = True
            return
        for line in self:
            if progress > i * alertInterval:
                print('Progress for %s: %s, temperature: %s, pressure: %s <%s%s>' % (layer.name, molecule.name, layer.T, layer.P, '*' * i, '-' * (20 - i)), end='\r')
                os.sys.stdout.flush()
                i += 1
            progress += 1
            xValues = np.arange(0, layer.distanceFromCenter, layer.resolution)
            hwRatio = line.lorentzHW / line.gaussianHW
            if hwRatio < .01:
                rightCurve = ls.gaussianLineShape(line.gaussianHW, xValues)
                trackGauss += 1
            elif hwRatio > 100:
                rightCurve = ls.lorentzLineShape(line.lorentzHW, xValues)
                trackLorentz += 1
            else:
                rightCurve = ls.pseudoVoigtShape(line.gaussianHW, line.lorentzHW, xValues)
                trackVoigt += 1
            intensity = pyradIntensity.intensityFactor(line.intensity, line.broadenedLine,
                                                       layer.T, line.lowerEnergy, self.q[layer.T], self.q296)
            arrayIndex = int((line.wavenumber - layer.rangeMin) / layer.resolution)
            arrayLength = len(crossSection) - 1
            if isBetween(arrayIndex, 0, arrayLength):
                crossSection[arrayIndex] = crossSection[arrayIndex] + rightCurve[0] * intensity
            for dx in range(1, len(rightCurve) - 1):
                rightIndex = arrayIndex + dx
                leftIndex = arrayIndex - dx
                if isBetween(rightIndex, 0, arrayLength):
                    crossSection[rightIndex] += rightCurve[dx] * intensity
                if isBetween(leftIndex, 0, arrayLength):
                    crossSection[leftIndex] += rightCurve[dx] * intensity
        self.crossSection = interpolateArray(self.xAxis,
                                             np.linspace(self.rangeMin, self.rangeMax,
                                                         (self.rangeMax - self.rangeMin) / self.resolution,
                                                         endpoint=True),
                                             crossSection)
        print('\ngaussian only: %s\t lorentz only: %s\t voigt: %s\n' % (trackGauss, trackLorentz, trackVoigt), end='\r')
        self.progressCrossSection = True

    def createLineSurvey(self):
        print('Creating line survey for %s' % self.name)
        molecule = self.molecule
        layer = molecule.layer
        progress = 0
        i = 1
        alertInterval = int(len(self) / 20)
        lineSurvey = np.zeros(int((self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION))
        for line in self:
            if progress > i * alertInterval:
                print('Progress for %s <%s%s>' % (molecule.name, '*' * i, '-' * (20 - i)), end='\r', flush=True)
                i += 1
            progress += 1
            intensity = line.intensity
            arrayIndex = int((line.wavenumber - layer.rangeMin) / layer.resolution)
            arrayLength = len(lineSurvey) - 1
            if isBetween(arrayIndex, 0, arrayLength):
                lineSurvey[arrayIndex] = lineSurvey[arrayIndex] + intensity
        self.lineSurvey = lineSurvey
        return self.lineSurvey

    def linelist(self):
        lines = []
        for line in self:
            lines.append(line)
        return lines

    def planck(self, temperature):
        return self.layer.planck(temperature)

    def transmission(self, surfaceSpectrum):
        transmitted = self.transmittance * surfaceSpectrum
        emitted = self.emittance * self.planck(self.T)
        return transmitted + emitted

    def returnCopy(self):
        return self.linelist()


class Molecule(list):
    def __init__(self, shortNameOrMolNum, layer, isotopeDepth=1, **abundance):
        super(Molecule, self).__init__(self)
        self.layer = layer
        self.crossSection = np.copy(layer.crossSection)
        self.isotopeDepth = isotopeDepth
        self.concText = ''
        self.concentration = 0
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
        self.progressCrossSection = False
        for key in abundance:
            if key == 'ppm':
                self.setPPM(abundance[key])
            elif key == 'ppb':
                self.setPPB(abundance[key])
            elif key == 'percentage' or key == 'perc' or key == '%':
                self.setPercentage(abundance[key])
            elif key == 'concentration':
                self.setConcentration(abundance[key])
            else:
                print('Invalid concentration type. Use ppm, ppb, percentage, or concentration.')

    def __str__(self):
        return '%s: %s' % (self.name, self.concText)

    def __bool__(self):
        return True

    def returnCopy(self):
        valueUnit = self.concText.split()
        tempDict = {valueUnit[1]: float(valueUnit[0])}
        newMolecule = Molecule(self.name, self.layer, isotopeDepth=int(self.isotopeDepth), **tempDict)
        newMolecule.getData(verbose=False)
        return newMolecule

    def setPercentage(self, percentage):
        self.concentration = percentage / 100
        self.concText = '%s %%' % percentage
        resetCrossSection(self)

    def setPPM(self, ppm):
        self.concentration = ppm * 10**-6
        self.concText = '%s ppm' % ppm
        resetCrossSection(self)

    def setPPB(self, ppb):
        self.concentration = ppb * 10**-8
        self.concText = '%s ppb' % ppb
        resetCrossSection(self)

    def setConcentration(self, concentration):
        self.setPPM(concentration * 1E6)
        resetCrossSection(self)

    def getData(self, verbose=True):
        for isotope in self:
            isotope.getData(verbose)

    def createCrossSection(self):
        tempAxis = np.zeros(int((self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION))
        for isotope in self:
            tempAxis += getCrossSection(isotope)
        self.progressCrossSection = True
        self.crossSection = tempAxis

    def planck(self, temperature):
        return self.layer.planck(temperature)

    def transmission(self, surfaceSpectrum):
        transmitted = self.transmittance * surfaceSpectrum
        emitted = self.emittance * self.planck(self.T)
        return transmitted + emitted

    @property
    def absCoef(self):
        return self.crossSection * self.concentration * self.layer.P / 1E4 / k / self.layer.T

    @property
    def transmittance(self):
        return np.exp(-self.absCoef * self.layer.depth)

    @property
    def lineSurvey(self):
        tempAxis = np.zeros(int((self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION))
        for isotope in self:
            tempAxis += isotope.lineSurvey
        return tempAxis

    @property
    def absorbance(self):
        return np.log10(1 / self.transmittance)

    @property
    def emissivity(self):
        return 1 - self.transmittance

    @property
    def emittance(self):
        return self.emissivity

    @property
    def P(self):
        return self.layer.P

    @property
    def yAxis(self):
        return np.copy(self.layer.yAxis)

    @property
    def xAxis(self):
        return np.copy(self.layer.xAxis)

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
    def molarMass(self):
        return self[0].molmass


class Layer(list):
    hasAtmosphere = False

    def __init__(self, depth, T, P, rangeMin, rangeMax, height=0.0, atmosphere=None, name='', dynamicResolution=True):
        super(Layer, self).__init__(self)
        self.rangeMin = rangeMin
        self.rangeMax = rangeMax
        self.T = T
        self.P = P
        self.depth = depth
        self.distanceFromCenter = self.P / 1013.25 * 5
        self.effectiveRangeMin = max(self.rangeMin - self.distanceFromCenter, 0)
        self.effectiveRangeMax = self.rangeMax + self.distanceFromCenter
        self.dynamicResolution = dynamicResolution
        self.surfaceSpectrum = None
        self.height = height
        if not dynamicResolution:
            self.resolution = utils.BASE_RESOLUTION
        else:
            self.resolution = max(10**int(np.log10((self.P / 1013.25))) * .01, utils.BASE_RESOLUTION)
        if not atmosphere:
            if not Layer.hasAtmosphere:
                self.atmosphere = Atmosphere('generic')
                Layer.hasAtmosphere = self.atmosphere
            else:
                self.atmosphere = Layer.hasAtmosphere
        else:
            self.atmosphere = atmosphere
            self.hasAtmosphere = atmosphere
        self.crossSection = np.zeros(int((rangeMax - rangeMin) / utils.BASE_RESOLUTION))
        self.progressCrossSection = False
        if not name:
            name = '%s' % self.atmosphere.nextLayerName()
        self.atmosphere.append(self)
        self.name = name

    def __str__(self):
        return '%s; %s' % (self.name, '; '.join(str(m) for m in self))

    def __bool__(self):
        return True

    def createCrossSection(self):
        tempAxis = np.zeros(int((self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION))
        for molecule in self:
            tempAxis += getCrossSection(molecule)
        self.progressCrossSection = True
        self.crossSection = tempAxis

    def returnMolecule(self, name):
        for m in self:
            if m.name == name:
                return m
        return False

    @property
    def meanHeight(self):
        return self.height + .5 * self.depth

    @property
    def lineSurvey(self):
        tempAxis = np.zeros(int((self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION))
        for molecule in self:
            tempAxis += molecule.lineSurvey
        return tempAxis

    @property
    def yAxis(self):
        return np.zeros(int((self.rangeMax - self.rangeMin) / self.resolution))

    @property
    def xAxis(self):
        return np.linspace(self.rangeMin, self.rangeMax, (self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION,
                           endpoint=True)

    @property
    def absCoef(self):
        tempAxis = np.zeros(int((self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION))
        for molecule in self:
            tempAxis += getAbsCoef(molecule)
        return tempAxis

    @property
    def transmittance(self):
        return np.exp(-self.absCoef * self.depth)

    @property
    def absorbance(self):
        return np.log10(1 / self.transmittance)

    @property
    def title(self):
        return '%s\nP: %smBars; T: %sK; depth: %scm' % (str(self), self.P, self.T, self.depth)

    @property
    def emissivity(self):
        return 1 - self.transmittance

    @property
    def emittance(self):
        return self.emissivity

    @property
    def molarMass(self):
        mass = 0
        for mol in self:
            mass += mol.molarMass * mol.concentration
        return mass

    @property
    def density(self):
        return self.P * 100 / self.T / self.specGasConstant

    @property
    def mass(self):
        return self.density * self.depth / 100

    @property
    def specGasConstant(self):
        return R * 1000 / self.molarMass

    @property
    def temperatureAtHeight(self):
        return int(self.atmosphere.temperatureAtHeight(self.meanHeight))

    @property
    def pressureAtHeight(self):
        return self.atmosphere.pressureAtHeight(self.meanHeight)

    @property
    def densityAtHeight(self):
        return self.atmosphere.densityAtHeight(self.temperatureAtHeight, self.specGasConstant)

    def changeRange(self, rangeMin, rangeMax):
        self.rangeMin = rangeMin
        self.rangeMax = rangeMax
        self.effectiveRangeMax = self.rangeMax + self.distanceFromCenter
        self.effectiveRangeMin = max(self.rangeMin - self.distanceFromCenter, 0)
        resetData(self)

    def changeTemperature(self, temperature):
        self.T = temperature
        resetCrossSection(self)

    def changePressure(self, pressure):
        self.P = pressure
        self.distanceFromCenter = self.P / 1013.25 * 5
        if not self.dynamicResolution:
            self.resolution = utils.BASE_RESOLUTION
        else:
            self.resolution = max(10**int(np.log10((self.P / 1013.25))) * .01, utils.BASE_RESOLUTION)
        resetData(self)

    def changeDepth(self, depth):
        self.depth = depth

    def addMolecule(self, name, isotopeDepth=1, **abundance):
        molecule = Molecule(name, self, isotopeDepth, **abundance)
        self.append(molecule)
        if totalConcentration(self) > 1:
            print('**Warning : Concentrations exceed 1.')
        molecule.getData()
        return molecule

    def returnCopy(self, name=None):
        if not name:
            name = self.atmosphere.nextLayerName()
        newCopy = Layer(self.depth, self.T, self.P, self.rangeMin, self.rangeMax, height=self.height,
                        atmosphere=self.atmosphere, name=name, dynamicResolution=self.dynamicResolution)
        for molecule in self:
            newCopy.addMolecule(molecule.name, molecule.isotopeDepth, concentration=molecule.concentration)
        print('%s copied to %s' % (self.name, newCopy.name))
        return newCopy

    def returnMoleculeObjects(self):
        moleculeList = []
        for m in self:
            moleculeList.append(m)
        return moleculeList

    def planck(self, temperature):
        return pyradPlanck.planckWavenumber(self.xAxis, temperature)

    def transmission(self, surfaceSpectrum):
        transmitted = self.transmittance * surfaceSpectrum
        emitted = self.emittance * self.planck(self.T)
        return transmitted + emitted


class Atmosphere(list):
    def __init__(self, name, rangeMin=0, rangeMax=0, planet=None):
        super().__init__(self)
        self.name = name
        self.rangeMin = rangeMin
        self.rangeMax = rangeMax
        self.planet = planet

    def __str__(self):
        return self.name

    def __bool__(self):
        return True

    def addLayer(self, depth, T, P, rangeMin, rangeMax, name=None, dynamicResolution=True, height=0.0):
        if not name:
            name = self.nextLayerName()
        newLayer = Layer(depth, T, P, rangeMin, rangeMax, atmosphere=self, name=name,
                         dynamicResolution=dynamicResolution, height=height)
        return newLayer

    def nextLayerName(self):
        return 'Layer %s' % (len(self) + 1)

    def returnLayerNames(self):
        tempList = []
        for layer in self:
            tempList.append(layer.name)
        return tempList

    def returnLayerObjects(self):
        if len(self) == 0:
            return False
        tempList = []
        for layer in self:
            tempList.append(layer)
        return tempList

    def temperatureAtHeight(self, height):
        ruleList = self.planet.returnApplicableRules(height, 'temperature')
        if len(ruleList) > 1:
            print('Multiple rules found for height %s: %s' % (height, ruleList))
        rule = ruleList[0]
        temperature = rule.rateFunction(rule.baseValue, rule.baseHeight, rule.rate, height)
        return temperature

    def pressureAtHeight(self, height):
        height = height / 100
        ruleList = self.planet.returnApplicableRules(height, 'temperature')
        if len(ruleList) > 1:
            print('Multiple rules found for height %s: %s' % (height, ruleList))
        rule = ruleList[0]
        temperature = self.temperatureAtHeight(height)
        if rule.rate != 0:
            pressure = rule.basePressure * (rule.baseValue / temperature) ** \
                  (self.planet.gravity * self.planet.molarMass / R / rule.rate)
        else:
            pressure = rule.basePressure * \
                       np.exp(-self.planet.gravity * self.planet.molarMass *
                              (height - rule.baseHeight) / R / rule.baseValue)
        return pressure

    def compositionAtHeight(self, height, molecule):
        ruleList = self.planet.returnApplicableRules(height, 'composition')
        if not ruleList:
            return molecule.concentration
        for rule in ruleList:
            if rule.molecule.name == molecule.name:
                concentration = rule.rateFunction(rule.baseValue, rule.baseHeight, rule.rate, height)
                return concentration
        return molecule.concentration

    def densityAtHeight(self, height, gasConstant):
        temperature = self.temperatureAtHeight(height)
        pressure = self.pressureAtHeight(height)
        density = pressure * 100 / gasConstant / temperature
        return density


class Planet:
    def __init__(self, name, pressure, temperature, molarMass=0.0289644,
                 gravity=9.80665, rangeMin=0, rangeMax=2000, initialThickness=100):
        self.name = name
        self.mass = 0
        self.gravity = gravity
        self.molarMass = molarMass
        self.radius = 0
        self.surfacePressure = pressure
        self.surfaceTemperature = temperature
        self.atmosphereRules = []
        self.rangeMin = rangeMin
        self.rangeMax = rangeMax
        self.atmosphere = Atmosphere("%s's atmosphere" % self.name, planet=self)
        self.initialLayer =  \
            self.atmosphere.addLayer(initialThickness, temperature, pressure, rangeMin, rangeMax, height=0)

    def setMass(self, mass):
        self.mass = mass

    def addLapseRate(self, name, baseHeight, baseTemperature, finalHeight, finalValue, basePressure, rateFunction=linear):
        self.atmosphereRules.append(AtmosphereRule(
            name, baseHeight, baseTemperature, finalHeight, finalValue, 'temperature',
            basePressure, rateFunction=rateFunction))

    def addCompositionRate(self, name, baseHeight, baseValue, finalHeight, rate, molecule, rateFunction=linear):
        self.atmosphereRules.append(AtmosphereRule(
            name, baseHeight, baseValue, finalHeight, rate, ruleType='composition',
            molecule=molecule, rateFunction=rateFunction))

    def setRadius(self, radius):
        self.radius = radius

    def setGravity(self, gravity):
        self.gravity = gravity

    def setSurfacePressure(self, pressure):
        self.surfacePressure = pressure

    def setSurfaceTemperature(self, temperature):
        self.surfaceTemperature = temperature

    def setSurfaceComposition(self):
        pass

    def returnApplicableRules(self, height, ruleType):
        ruleList = []
        for rule in self.atmosphereRules:
            if rule.isInRange(height, ruleType):
                ruleList.append(rule)
        if not ruleList:
            return False
        return ruleList

    def temperatureAtHeight(self, height):
        return self.atmosphere.temperatureAtHeight(height)

    def pressureAtHeight(self, height):
        return self.atmosphere.pressureAtHeight(height)

    def densityAtHeight(self, height):
        return self.atmosphere.densityAtHeight(height)

    def compositionAtHeight(self, height, molecule):
        return self.atmosphere.compositionAtHeight(height, molecule)

    def addMolecule(self, name, isotopeDepth=1, **abundance):
        molecule = self.initialLayer.addMolecule(name, isotopeDepth, **abundance)
        return molecule

    def sliceAtm(self):
        acceptSetup = False
        height = 0
        while not acceptSetup:
            mass = self.initialLayer.mass
            layer = self.initialLayer.returnCopy(name='temp layer')
            print('p %s, T %s, depth %s, height %s' % (layer.P, layer.T, layer.depth, layer.height))
            heightList = [height]
            tempList = [layer.T]
            pressureList = [layer.P]
            depthList = [layer.depth]
            changePoints = self.atmChangePoints
            while layer.height + layer.depth < 9000000:
                print('temp: %s, press: %s, baseHeight: %skm, meanHeight: %s, depth: %s, tarMass: %s, layMass: %s' % (
                       int(layer.T), layer.P, layer.height / 100000, layer.meanHeight, layer.depth, mass, layer.mass))
                newHeight = heightList[-1] + depthList[-1]
                heightList.append(newHeight)
                tempList.append(self.temperatureAtHeight(newHeight))
                pressureList.append(self.pressureAtHeight(newHeight))
                layer.P = pressureList[-1]
                layer.T = tempList[-1]
                layer.height = heightList[-1]
                layer.depth = (mass / layer.density * 100)
                '''if layer.height + layer.depth > changePoints[0]:
                    tempHeight = changePoints.pop(0)
                    layer.depth = tempHeight - layer.height
                    while changePoints:
                        if changePoints[0] == tempHeight:
                            changePoints.pop(0)
                        else:
                            break'''
                if layer.depth != 0:
                    depthList.append(layer.depth)
            print('Total # of layers: %s' % len(heightList))
            isValid = False
            while not isValid:
                self.atmosphere.pop()
                userInput = input('Accept the current slicing (y/n): ')
                if userInput.lower()[0] == 'n':
                    isValid = True
                    validNumber = False
                    while not validNumber:
                        print('Enter the new initial depth. Larger depth will decrease number of layers, smaller will increase')
                        userNumber = input('Current depth is %scm:' % utils.limeText(self.initialLayer.depth))
                        try:
                            newDepth = float(userNumber)
                            self.initialLayer.depth = newDepth
                            validNumber = True
                        except ValueError:
                            print('Invalid number.')
                elif userInput.lower()[0] == 'y':
                    isValid = True
                    acceptSetup = True
        # print(len(self.atmosphere), self.atmosphere[0].name)
        heightList.pop(0)
        depthList.pop(0)
        for i in range(0, len(heightList) - 1):
            height = heightList[i]
            depth = depthList[i]
            newLayer = self.initialLayer.returnCopy()
            newLayer.height = height
            newLayer.depth = depth
            newLayer.T = int(self.temperatureAtHeight(newLayer.meanHeight))
            newLayer.P = self.pressureAtHeight(newLayer.meanHeight)
            for molecule in newLayer:
                molecule.concentration = self.compositionAtHeight(newLayer.meanHeight, molecule)
        surfaceSpectrum = pyradPlanck.planckWavenumber(self.initialLayer.xAxis, self.surfaceTemperature)

        for layer in self.atmosphere:
            layer.createCrossSection()
            transmission = layer.transmission(surfaceSpectrum)
            surfaceSpectrum = transmission
        plt.plot(self.initialLayer.xAxis, smooth(surfaceSpectrum, 100), linewidth=.5)
        plt.plot(self.initialLayer.xAxis, pyradPlanck.planckWavenumber(self.initialLayer.xAxis, self.surfaceTemperature), linewidth=.5)
        plt.plot(self.initialLayer.xAxis,
                 pyradPlanck.planckWavenumber(self.initialLayer.xAxis, 280), linewidth=.5)
        plt.plot(self.initialLayer.xAxis,
                 pyradPlanck.planckWavenumber(self.initialLayer.xAxis, 260), linewidth=.5)
        plt.plot(self.initialLayer.xAxis,
                 pyradPlanck.planckWavenumber(self.initialLayer.xAxis, 240), linewidth=.5)
        plt.show()

    @property
    def yAxis(self):
        return np.zeros(int((self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION))

    @property
    def xAxis(self):
        return np.linspace(self.rangeMin, self.rangeMax, (self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION,
                           endpoint=True)

    @property
    def atmChangePoints(self):
        changePointList = []
        for rule in self.atmosphereRules:
            changePointList.append(rule.finalHeight)
        return sorted(changePointList)


class AtmosphereRule:
    def __init__(self, name, baseHeight, baseValue, finalHeight, finalValue, ruleType, basePressure=0,
                 molecule=None, rateFunction=linear):
        self.name = name
        self.baseHeight = baseHeight
        self.baseValue = baseValue
        self.finalHeight = finalHeight
        self.finalValue = finalValue
        if rateFunction == linear:
            self.rate = (self.finalValue - self.baseValue) / (self.finalHeight - self.baseHeight)
        self.ruleType = ruleType
        self.molecule = molecule
        self.basePressure = basePressure
        self.rateFunction = rateFunction

    def isInRange(self, height, ruleType):
        return self.baseHeight < height <= self.finalHeight and self.ruleType == ruleType


def returnPlot(obj, propertyToPlot):
    if propertyToPlot == "transmittance":
        yAxis = getTransmittance(obj), 1
    elif propertyToPlot == 'absorption coefficient':
        yAxis = getAbsCoef(obj), 0
    elif propertyToPlot == 'cross section':
        yAxis = getCrossSection(obj), 0
    elif propertyToPlot == 'absorbance':
        yAxis = getAbsorbance(obj), 0
    elif propertyToPlot == 'optical depth':
        yAxis = getOpticalDepth(obj), 0
    elif propertyToPlot == 'line survey':
        yAxis = obj.lineSurvey, 0
    else:
        return False
    return yAxis


def isBetween(test, minValue, maxValue):
    if test >= minValue:
        if test <= maxValue:
            return True
    return False


def plot(propertyToPlot, title, plotList, fill=False):
    plt.figure(figsize=(10, 6), dpi=80)
    plt.subplot(111, facecolor='xkcd:dark grey')
    plt.xlabel('wavenumber cm-1')
    plt.margins(0.01)
    plt.subplots_adjust(left=.07, bottom=.08, right=.97, top=.90)
    plt.ylabel(propertyToPlot)
    if propertyToPlot == 'line survey':
        plt.yscale('log')
    plt.grid('grey', linewidth=.5, linestyle=':')
    plt.title('%s' % title)
    handles = []
    linewidth = 1.2
    alpha = .7
    for singlePlot, color in zip(plotList, COLOR_LIST):
        yAxis, fillAxis = returnPlot(singlePlot, propertyToPlot)
        fig, = plt.plot(singlePlot.xAxis,  smooth(yAxis, 50), linewidth=linewidth, alpha=alpha, color=color,
                        label='%s' % singlePlot.name)
        handles.append(fig)
        plt.fill_between(singlePlot.xAxis, fillAxis, yAxis, color=color, alpha=.3 * fill)
        linewidth = .7
        alpha = .5
    legend = plt.legend(handles=handles, frameon=False)
    text = legend.get_texts()
    plt.setp(text, color='w')
    plt.show()


def plotSpectrum(layer=None, title=None, rangeMin=None, rangeMax=None, objList=None, surfaceSpectrum=None,
                 planckTemperatureList=None, planckType='wavenumber'):
    plt.figure(figsize=(10, 6), dpi=80)
    plt.subplot(111, facecolor='xkcd:dark grey')
    plt.margins(0.01)
    plt.subplots_adjust(left=.07, bottom=.08, right=.97, top=.90)
    if layer:
        rangeMin = layer.rangeMin
        rangeMax = layer.rangeMax
        title = layer.title
    if planckType == 'wavenumber':
        plt.xlabel('wavenumber cm-1')
        plt.ylabel('Radiance Wm-2sr-1(cm-1)-1')
        planckFunction = pyradPlanck.planckWavenumber
        xAxis = np.linspace(rangeMin, rangeMax, (rangeMax - rangeMin) / utils.BASE_RESOLUTION)
    elif planckType == 'Hz':
        plt.xlabel('Hertz')
        plt.ylabel('Radiance Wm-2sr-1Hz-1')
        planckFunction = pyradPlanck.planckHz
        xAxis = np.linspace(rangeMin, rangeMax, 1000)
    else:
        plt.xlabel('wavelength um')
        plt.ylabel('Radiance Wm-2sr-1um-1')
        planckFunction = pyradPlanck.planckWavelength
        xAxis = np.linspace(rangeMin, rangeMax, (rangeMax - rangeMin) / utils.BASE_RESOLUTION)
    plt.title('%s' % title)
    handles = []
    blue = .3
    red = 1
    green = .6
    dr = -.15
    db = .15
    dg = .15
    if not rangeMax:
        xAxis = layer.xAxis
    for temperature in planckTemperatureList:
        yAxis = planckFunction(xAxis, float(temperature))
        fig, = plt.plot(xAxis, yAxis, linewidth=.75, color=(red, green, blue),
                        linestyle=':', label='%sK : %sWm-2' %
                        (temperature, round(integrateSpectrum(yAxis, res=(rangeMax - rangeMin) / len(yAxis)), 2)))
        handles.append(fig)
        if red + dr < 0 or red + dr > 1:
            dr *= -1
        if green + dg > 1 or green + dg < 0:
            dg *= -1
        if blue + db > 1 or blue + db < 0:
            db *= -1
        red += dr
        green += dg
        blue += db
        if red < .3 and green < .3 and blue < .3:
            green += .5
            blue += .2
        if red < .3 and green < .3:
            green += .4
    if objList:
        alpha = .7
        linewidth = 1.2
        for obj, color in zip(objList, COLOR_LIST):
            yAxis = obj.transmission(surfaceSpectrum)
            fig, = plt.plot(layer.xAxis, yAxis, linewidth=linewidth,
                            alpha=alpha, color=color,
                            label='%s : %sWm-2' % (obj.name, round(integrateSpectrum(yAxis, pi), 2)))
            handles.append(fig)
            alpha = .5
            linewidth = 1
    legend = plt.legend(handles=handles, frameon=False)
    text = legend.get_texts()
    plt.setp(text, color='w')
    plt.show()


def cacheCurves():
    ls.writeCacheToFile()


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
                     49: {1: 124, 2: 125},
                     901: {1: 901}}

COLOR_LIST = ['xkcd:white',
              'xkcd:bright orange',
              'xkcd:seafoam green',
              'xkcd:bright blue',
              'xkcd:salmon',
              'xkcd:light violet',
              'xkcd:green yellow']

VERSION = utils.VERSION

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
               'cs': 46, 'so3': 47, 'c2n2': 48, 'cocl2': 49,
               'ar': 901}


if __name__ == 'main':
    pyradInteractive.menuMain()
