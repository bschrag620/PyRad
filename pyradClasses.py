import pyradUtilities as utils
import pyradLineshape as ls
import pyradIntensity
import pyradPlanck
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import time
import gc


c = 299792458.0
k = 1.38064852E-23
h = 6.62607004e-34
pi = 3.141592653589793
R = 8.3144598
t0 = 296
p0 = 1013.25
avo = 6.022140857E23
sb = 5.67E-8


settings = utils.Settings('high')
theme = utils.Theme()


def stefanB(power):
    return (power / sb) ** .25


def reduceRes(array, finalRes=1.0):
    n = int(finalRes / settings.baseRes)
    length = int(len(array) * settings.baseRes / finalRes)
    newArray = np.zeros(length)
    for m in range(0, length):
        for i in range(0, n):
            newArray[m] += array[m * n + i] / n
    return newArray


def strToBool(v):
    return 'true' in v.lower()


def returnHMS(time):
    hours = int(time / 3600)
    minutes = int((time - hours * 3600) / 60)
    seconds = time - hours * 3600 - minutes * 60
    return hours, minutes, seconds


def loadEmptyPlanet(folderPath, planet=None, verify=False):
    values = utils.readCompleteProfile(folderPath)
    if planet is None:
        planet = Planet(values['name'], float(values['surfacePressure']), int(values['surfaceTemperature']), float(values['maxHeight']),
                        rangeMin=int(values['rangeMin']), rangeMax=int(values['rangeMax']), initialThickness=int(float(values['initialDepth'])),
                        gravity=float(values['gravity']), res=int(float(values['res'])))
        for mol in values['molList'].split(','):
            planet.moleculeList.append(mol)
    if not utils.profileComplete(folderPath):
        planet.processLayers(verify=verify)
    fileLength = utils.profileLength(folderPath)
    while len(planet.atmosphere) > 0:
        planet.atmosphere.pop()
    for i in range(1, fileLength + 1):
        lP = utils.readPlanetProfile(folderPath, i, fileLength)
        layer = Layer(lP['depth'], lP['T'], lP['P'], lP['rangeMin'], lP['rangeMax'], height=lP['height'],
                      name=lP['name'])
        layer.absorptionCoefficient = np.asarray(lP['absCoef'])
        layer.progressAbsCoef = True
        planet.atmosphere.append(layer)
    print('')
    planet.progressProfileLoaded = True
    return planet


def readTransmissionFromFile(requestedHeight, folderPath, direction):
    completedValues = utils.readCompleteTransmission(folderPath)
    heightList = completedValues['heightList']
    maxHeight = float(completedValues['maxHeight']) * 100000
    i = 0

    if direction == 'down':
        height = 0
        while height < requestedHeight and i < len(heightList) - 1:
            height = float(heightList[i])
            i += 1
    else:
        height = maxHeight
        heightList.reverse()
        requestedHeight = min(requestedHeight, maxHeight)
        while height > requestedHeight and i < len(heightList) - 1:
            height = float(heightList[i])
            i += 1
    targetIndex = heightList.index(height)
    fileName = 'trans looking %s-%s.pyr' % (direction, targetIndex)
    print('retreiving transmission from: %s' % fileName)
    transmissionValues = utils.readTransmissionValues(fileName, folderPath)
    for key in transmissionValues:
        transmissionValues[key] = np.asarray(transmissionValues[key])
    transmissionValues['rangeMin'] = float(completedValues['rangeMin'])
    transmissionValues['rangeMax'] = float(completedValues['rangeMax'])
    transmissionValues['surfaceTemperature'] = float(completedValues['surfaceTemperature'])
    transmissionValues['molList'] = completedValues['molList']
    transmissionValues['surfaceEffEmissivity'] = float(completedValues['surfaceEffEmissivity']),
    transmissionValues['res'] = float(completedValues['res'])
    transmissionValues['surfacePower'] = float(completedValues['surfacePower'])
    return transmissionValues


def processTransmissionBySingleLayer(folderPath, res=1):

    values = utils.readCompleteProfile(folderPath)
    planet = Planet(values['name'], float(values['surfacePressure']), int(float(values['surfaceTemperature'])),
                    float(values['maxHeight']),
                    rangeMin=int(values['rangeMin']), rangeMax=int(values['rangeMax']),
                    initialThickness=int(float(values['initialDepth'])),
                    gravity=float(values['gravity']), res=int(float(values['res'])), surfaceEffEmissivity=float(values['surfaceEffEmissivity']))
    for mol in values['molList'].split(','):
        planet.moleculeList.append(mol)

    xAxis = None
    fileLength = utils.profileLength(folderPath)
    heightList = [0]
    # making a generic layer and molecules
    layer = Layer(0, 0, 0, 0, 0)

    spectrumDict = {}
    inputDict = {}
    
    for i in range(1, fileLength + 1):
        layerProfile = utils.readPlanetProfile(folderPath, i, fileLength)
        layer.depth = layerProfile['depth']
        layer.P = layerProfile['P']
        layer.T = layerProfile['T']
        layer.rangeMin = layerProfile['rangeMin']
        layer.rangeMax = layerProfile['rangeMax']
        layer.height = layerProfile['height']
        layer.name = layerProfile['name']
        layer.absorptionCoefficient = np.asarray(layerProfile['layer absCoef'])
        layer.progressAbsCoef = True

        if xAxis is None:
            xAxis = np.linspace(planet.rangeMin, planet.rangeMax, len(layer.absorptionCoefficient))
            spectrumDict['temperature'] = planet.surfaceTemperature,
            spectrumDict['pressure'] = planet.surfacePressure,
            spectrumDict['depth'] = 'surface',
            spectrumDict['meanHeight'] = 0
            inputSpectrum = pyradPlanck.planckWavenumber(xAxis, planet.surfaceTemperature) * planet.effEmissivity
            surfacePower = integrateSpectrum(inputSpectrum, pi,.01)
            spectrumDict['layer transmission'] = reduceRes(inputSpectrum, finalRes=res)
            spectrumDict['surfacePower'] = surfacePower
            inputDict['layer'] = inputSpectrum
            for molecule in planet.moleculeList:
                spectrumDict['%s transmission' % molecule] = reduceRes(inputSpectrum, finalRes=res)
                spectrumDict['%s concentration' % molecule] = layerProfile['%s concentration' % molecule]
                spectrumDict['%s power' % molecule] = integrateSpectrum(inputSpectrum, pi, .01)
                inputDict[molecule] = inputSpectrum
            utils.writePlanetTransmission(folderPath, 0, spectrumDict, 'down', 0)
        spectrumDict['temperature'] = layer.T,
        spectrumDict['pressure'] = layer.P,
        spectrumDict['depth'] = layer.depth,
        spectrumDict['meanHeight'] = layer.meanHeight
        inputSpectrum = inputDict['layer']
        transmittedSpectrum = layer.transmission(inputSpectrum)
        spectrumDict.update({'layer transmission': reduceRes(transmittedSpectrum, finalRes=res),
                             'layer effective emissivity': layer.effectiveEmissivity,
                             'layer power': integrateSpectrum(transmittedSpectrum, pi, .01)})
        inputDict['layer'] = transmittedSpectrum

        for molecule in planet.moleculeList:
            inputSpectrum = inputDict[molecule]
            layer.absorptionCoefficient = np.asarray(layerProfile['%s absCoef' % molecule])
            transmittedSpectrum = layer.transmission(inputSpectrum)
            spectrumDict.update({'%s transmission' % molecule: reduceRes(transmittedSpectrum),
                                 '%s effective emissivity' % molecule: layer.effectiveEmissivity,
                                 '%s power' % molecule: integrateSpectrum(transmittedSpectrum, pi, .01),
                                 '%s concentration' % molecule: layerProfile['%s concentration' % molecule]})
            inputDict[molecule] = transmittedSpectrum

        utils.writePlanetTransmission(folderPath, layer.meanHeight, spectrumDict, 'down', i)
        heightList.append(layer.meanHeight)
        gc.collect()

    heightList.append(planet.maxHeight)

    # with transmission from surface upward processed, do the same in reverse to get the transmission toward the surface
    # initial spectrum will be 2.7K for CMB
    inputDict = {}

    inputSpectrum = pyradPlanck.planckWavenumber(xAxis, 2.7)
    spectrumDict = {'layer transmission': reduceRes(inputSpectrum, finalRes=res)}
    surfacePower = integrateSpectrum(inputSpectrum, pi, .01)
    spectrumDict['surfacePower'] = surfacePower
    inputDict['layer'] = inputSpectrum
    for molecule in planet.moleculeList:
        spectrumDict['%s transmission' % molecule] = reduceRes(inputSpectrum, finalRes=res)
        inputDict[molecule] = inputSpectrum
    utils.writePlanetTransmission(folderPath, planet.maxHeight, spectrumDict, 'up', 0)

    for i in range(1, fileLength + 1):
        fileNumber = fileLength + 1 - i
        layerProfile = utils.readPlanetProfile(folderPath, fileNumber, fileLength)
        layer.depth = layerProfile['depth']
        layer.P = layerProfile['P']
        layer.T = layerProfile['T']
        layer.rangeMin = layerProfile['rangeMin']
        layer.rangeMax = layerProfile['rangeMax']
        layer.height = layerProfile['height']
        layer.name = layerProfile['name']
        layer.absorptionCoefficient = np.asarray(layerProfile['layer absCoef'])
        layer.progressAbsCoef = True
        inputSpectrum = inputDict['layer']
        transmittedSpectrum = layer.transmission(inputSpectrum)
        spectrumDict.update({'layer transmission': reduceRes(transmittedSpectrum, finalRes=res),
                        'layer effective emissivity': layer.effectiveEmissivity,
                        'layer power': integrateSpectrum(transmittedSpectrum, pi, .01)})
        inputDict['layer'] = transmittedSpectrum

        for molecule in planet.moleculeList:
            inputSpectrum = inputDict[molecule]
            layer.absorptionCoefficient = np.asarray(layerProfile['%s absCoef' % molecule])
            transmittedSpectrum = layer.transmission(inputSpectrum)
            spectrumDict.update({'%s transmission' % molecule: reduceRes(transmittedSpectrum),
                                 '%s effective emissivity' % molecule: layer.effectiveEmissivity,
                                 '%s power' % molecule: integrateSpectrum(transmittedSpectrum, pi, .01),
                                 '%s concentration' % molecule: layerProfile['%s concentration' % molecule]})
            inputDict[molecule] = transmittedSpectrum

        utils.writePlanetTransmission(folderPath, layer.meanHeight, spectrumDict, 'up', i)
        gc.collect()

    utils.profileWriteTransmissionComplete(folderPath, heightList)
    return


def createCustomPlanet(name):

    initialParameters = utils.parseCustomProfile(name)
    initialValues = initialParameters['initialValues']
    moleculeList = initialParameters['molecules']
    temperatureRuleList = initialParameters['temperatureRules']
    compositionRuleList = initialParameters['compositionRules']

    planet = Planet(initialParameters['name'], initialValues['surfacePressure'], initialValues['surfaceTemperature'],
                    initialValues['maxHeight'], rangeMin=initialValues['rangeMin'],
                    rangeMax=initialValues['rangeMax'], initialThickness=initialValues['initialDepth'],
                    gravity=initialValues['gravity'])

    for molecule in moleculeList:
        planet.addMolecule(molecule['name'], concentration=molecule['concentration'])

    for rule in temperatureRuleList:
        planet.addLapseRate(rule['name'], rule['finalHeight'], rule['finalValue'])

    for rule in compositionRuleList:
        initalLayer = planet.initialLayer
        molecule = initalLayer.returnMolecule(rule['moleculeName'])
        planet.addCompositionRate(rule['name'], rule['finalHeight'], rule['finalValue'], molecule)
    return planet


def yesOrNo(promptText):
    validInput = False
    while not validInput:
        userSelection = input(promptText)
        if not userSelection:
            pass
        elif userSelection.lower()[0] == 'y':
            return True
        elif userSelection.lower()[0] == 'n':
            return False
        else:
            print('Invalid choice.')


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


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
        if number in INERT_MOL_DATA:
            self.globalIsoNumber = MOLECULE_ID[molecule.name]
            self.shortName = molecule.name
            self.name = 'Inert molecule %s' % molecule.name
            self.molmass = INERT_MOL_DATA[self.globalIsoNumber]['mass']
        else:
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
        self.progressCrossSection = False

    def __del__(self):
        for line in self:
            line = None


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

    def getData(self, verbose=True):
        if 'Inert' in self.name:
            print('Inert molecule, no data.')
            return
        if verbose:
            print('Getting data for %s, %s, isotope %s' % (self.layer.name, self.molecule.name, self.globalIsoNumber), end='\r', flush=True)
        lineDict = utils.gatherData(self.globalIsoNumber, self.layer.effectiveRangeMin, self.layer.effectiveRangeMax)
        self.q = utils.getQData(self.globalIsoNumber)
        for line in lineDict:
            if lineDict[line]['intensity'] > settings.lineIntensityCutoff:
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
        crossSection = np.copy(self.yAxis)
        trackGauss = 0
        trackLorentz = 0
        trackVoigt = 0
        if len(self) == 0:
            self.progressCrossSection = True
            return
        for line in self:
            if progress > i * alertInterval and len(self) > 50:
                text = 'Progress %s: %s: <%s%s>' % (layer.name, molecule.name, '*' * i, '-' * (20 - i))
                print(text, end='\r', flush=True)
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
            if intensity > settings.lineIntensityCutoff:
                arrayIndex = int((line.wavenumber - layer.rangeMin) / layer.resolution)
                arrayLength = len(lineSurvey) - 1
                if isBetween(arrayIndex, 0, arrayLength):
                    lineSurvey[arrayIndex] = lineSurvey[arrayIndex] + intensity
        return lineSurvey

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

    @property
    def lineSurvey(self):
        return self.createLineSurvey()


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

    def __del__(self):
        for iso in self:
            iso = None

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
    def effectiveEmissivity(self):
        return np.average(self.emissivity)

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
            self.resolution = max(self.P / 1013.25 * .01, utils.BASE_RESOLUTION)
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
        self.absorptionCoefficient = np.zeros(int((rangeMax - rangeMin) / utils.BASE_RESOLUTION))
        self.progressCrossSection = False
        self.progressAbsCoef = False
        if not name:
            name = '%s' % self.atmosphere.nextLayerName()
        self.atmosphere.append(self)
        self.name = name

    def __str__(self):
        return '%s; %s' % (self.name, '; '.join(str(m) for m in self))

    def __bool__(self):
        return True

    def __del__(self):
        for molecule in self:
            molecule.__del__()

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
        return np.linspace(self.rangeMin, self.rangeMax, len(self.absorptionCoefficient),
                           endpoint=True)

    @property
    def absCoef(self):
        if not self.progressAbsCoef:
            for molecule in self:
                self.absorptionCoefficient += getAbsCoef(molecule)
        self.progressAbsCoef = True
        return self.absorptionCoefficient

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
    def effectiveEmissivity(self):
        return np.average(self.emissivity)

    @property
    def normalizedEmissivity(self):
        normalized = np.array([])
        emissivity = self.emissivity
        nonzeroIndex = np.nonzero(emissivity)
        if np.size(nonzeroIndex) == 0:
            return 0
        for i in nonzeroIndex:
            normalized = np.append(normalized, emissivity[i])
        return np.average(normalized)

    @property
    def molarMass(self):
        mass = 0
        for mol in self:
            if mol.name == 'h2o':
                pass
            else:
                mass += mol.molarMass * mol.concentration
        return mass

    @property
    def density(self):
        d = self.P * 100 / self.T / self.specGasConstant
        return d

    @property
    def mass(self):
        return self.density * self.depth / 100

    @property
    def specGasConstant(self):
        con = R * 1000 / self.molarMass
        return con

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
            self.resolution = max(self.P / 1013.25 * .01, utils.BASE_RESOLUTION)
        resetData(self)

    def changeDepth(self, depth):
        self.depth = depth

    def addMolecule(self, name, isotopeDepth=1, **abundance):
        molecule = Molecule(name, self, isotopeDepth, **abundance)
        self.append(molecule)
        if totalConcentration(self) > 1:
            print('**Warning : Concentrations exceed 1: total=%s' % totalConcentration(self))
            for molecule in self:
                print('concentration %s: %s' % (molecule.name, molecule.concentration))
            exit(1)
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

    def returnMoleculeNameList(self):
        nameList = []
        for m in self:
            nameList.append(m.name)
        return nameList

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

    def __del__(self):
        for layer in self:
            layer.__del__()

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
        ruleList = self.planet.returnApplicableRules(height, 'temperature')
        if len(ruleList) > 1:
            print('Multiple rules found for height %s: %s' % (height, ruleList))
        rule = ruleList[0]
        temperature = self.temperatureAtHeight(height)
        pressurePa = rule.basePressure * 100
        rateInMeters = rule.rate * 100
        if rule.rate != 0:
            pressure = pressurePa * (rule.baseValue / temperature) ** \
                  (self.planet.gravity * self.planet.molarMass / R / rateInMeters)
        else:
            changeInHeightMeters = (height - rule.baseHeight) / 100
            pressure = pressurePa * np.exp(-self.planet.gravity * self.planet.molarMass * changeInHeightMeters / R / rule.baseValue)
        return pressure / 100

    def compositionAtHeight(self, height, molecule):
        ruleList = self.planet.returnApplicableRules(height, 'composition')
        if not ruleList:
            return molecule.concentration
        for rule in ruleList:
            if rule.molecule.name == molecule.name:
                concentration = rule.rateFunction(rule.baseValue, rule.baseHeight, rule.rate, height)
                return concentration
        return molecule.concentration

    def densityAtHeight(self, height):
        gasConstant = self.planet.specGasConstant
        temperature = self.temperatureAtHeight(height)
        pressure = self.pressureAtHeight(height)
        density = pressure * 100 / gasConstant / temperature
        return density

    def moleculesAtHeight(self, height):
        density = self.densityAtHeight(height)
        return density / self.planet.molarMass / 1000 * avo


class Planet:
    def __init__(self, name, pressure, temperature, maxHeight, surfaceEffEmissivity=.971,
                 gravity=9.80665, rangeMin=0, rangeMax=2000, initialThickness=100, res=1):
        self.name = name
        self.gravity = gravity
        self.maxHeight = maxHeight * 100000
        self.effEmissivity = surfaceEffEmissivity
        self.surfacePressure = pressure
        self.surfaceTemperature = temperature
        self.atmosphereRules = []
        self.heightList = []
        self.depthList = []
        self.moleculeList = []
        self.rangeMin = rangeMin
        self.rangeMax = rangeMax
        self.setting = settings.setting
        self.atmosphere = Atmosphere("%s's atmosphere" % self.name, planet=self)
        self.res = res
        self.initialLayer =  \
            self.atmosphere.addLayer(initialThickness * 100, temperature, pressure, rangeMin, rangeMax,
                                     name='initial layer', height=0)
        self.progressProfileLoaded = False

    def __del__(self):
        self.atmosphere.__del__()

    def clearData(self):
        self.__del__()

    def addLapseRate(self, name, finalHeight, finalValue, rateFunction=linear):
        self.atmosphereRules.append(AtmosphereRule(name, finalHeight * 100000,
                                                   finalValue, 'temperature', self, rateFunction=rateFunction))

    def addCompositionRate(self, name, finalHeight, finalValue, molecule, rateFunction=linear):
        self.atmosphereRules.append(AtmosphereRule(name, finalHeight * 100000, finalValue, 'composition', self,
            molecule=molecule, rateFunction=rateFunction))

    def returnApplicableRules(self, height, ruleType):
        ruleList = []
        for rule in self.atmosphereRules:
            if rule.isInRange(height, ruleType):
                ruleList.append(rule)
        if not ruleList:
            return False
        return ruleList

    def returnAllRulesOfType(self, ruleType):
        ruleList = []
        for rule in self.atmosphereRules:
            if rule.ruleType == ruleType:
                ruleList.append(rule)
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
        self.moleculeList = self.initialLayer.returnMoleculeNameList()
        return molecule

    def sliceAtm(self, verify=True):
        acceptSetup = False
        layer = self.initialLayer
        initialTemp = layer.T
        initialPressure = layer.P
        initialHeight = layer.height
        initialDepth = layer.depth
        while not acceptSetup:
            mass = self.initialLayer.mass
            print('%s, p %s, T %s, depth %s, height %s' % (layer.name, initialPressure, initialTemp, initialDepth, initialHeight))
            self.heightList = [layer.height]
            tempList = [layer.T]
            pressureList = [layer.P]
            self.depthList = [layer.depth]
            while layer.height + layer.depth < self.maxHeight:
                print('%s: K: %s, mBar: %s, height: %skm, depth: %skm' % (
                       'Layer %s' % len(self.heightList), int(layer.T), round(layer.P, 2), round(layer.height / 100000, 2), round(layer.depth / 100000, 2)))
                newHeight = self.heightList[-1] + self.depthList[-1]
                tempList.append(self.temperatureAtHeight(newHeight))
                pressureList.append(self.pressureAtHeight(newHeight))
                layer.P = pressureList[-1]
                layer.T = tempList[-1]
                newDepth = (mass / layer.density * 100)
                if newHeight + newDepth > self.maxHeight:
                    newDepth = self.maxHeight - newHeight
                    self.heightList.append(newHeight)
                    self.depthList.append(newDepth)
                    print('%s: K: %s, mBar: %s, height: %skm, depth: %skm' % (
                        'Layer %s' % len(self.heightList), int(layer.T), round(layer.P, 2),
                        round(newHeight / 100000, 2), round(newDepth / 100000, 2)))
                    break
                else:
                    self.heightList.append(newHeight)
                    self.depthList.append(newDepth)
                    layer.height = newHeight
                    layer.depth = newDepth
            print('Total # of layers: %s' % len(self.heightList))
            print('Total # of absorption lines per layer: %s' % len(totalLineList(layer)))
            if verify:
                if yesOrNo('Accept the current slicing (y/n): '):
                    acceptSetup = True
                else:
                    validNumber = False
                    while not validNumber:
                        print('Enter the new initial depth. Larger depth will decrease number of layers, smaller will increase')
                        userNumber = input('Current depth is %s:' % utils.limeText('%sm' % (int(self.depthList[0]) / 100)))
                        try:
                            newDepth = float(userNumber)
                            self.initialLayer.depth = newDepth * 100
                            self.initialLayer.P = initialPressure
                            self.initialLayer.T = initialTemp
                            self.initialLayer.height = initialHeight
                            validNumber = True
                        except ValueError:
                            print('Invalid number.')
            else:
                acceptSetup = True
        return

    def processLayers(self, verify=True, moleculeSpecific=False, res=1):
        if utils.profileProgress(self.folderPath):
            self.sliceAtm(verify=False)
        if not self.heightList:
            self.sliceAtm(verify=verify)
        layer = self.initialLayer
        startPoint, timeStart = utils.profileProgress(self.folderPath)
        i = startPoint + 1
        totalProcessTime = timeStart
        for height, depth in zip(self.heightList[startPoint:], self.depthList[startPoint:]):
            layerProcessTimeStart = time.time()
            layer.name = 'layer %s_%s' % (i, len(self.heightList))
            layer.height = height
            layer.depth = depth
            layer.T = int(self.temperatureAtHeight(layer.meanHeight))
            layer.P = self.pressureAtHeight(layer.meanHeight)
            for molecule in layer:
                molecule.concentration = self.compositionAtHeight(layer.meanHeight, molecule)
            layer.createCrossSection()
            processTime = time.time() - layerProcessTimeStart
            utils.writePlanetProfile(self.folderPath, layer, processTime, self.moleculeList, moleculeSpecific=moleculeSpecific)
            totalProcessTime += processTime
            utils.profileWriteProgress(self.folderPath, i, len(self.heightList), totalProcessTime,
                                       moleculeSpecific, self.moleculeList, res)
            resetCrossSection(layer)
            layer.absorptionCoefficient = np.zeros(len(layer.crossSection))
            layer.progressAbsCoef = False
            i += 1
        utils.profileWriteComplete(self, i - 1, len(self.heightList), totalProcessTime, res, moleculeSpecific=moleculeSpecific)
        return

    def loadProfile(self, verify=True, moleculeSpecific=False):
        if not utils.profileComplete(self.folderPath):
            self.processLayers(verify=verify, moleculeSpecific=moleculeSpecific)
        loadEmptyPlanet(self.folderPath, self, verify=verify)

    def processTransmission(self, height, direction='down', verify=True, moleculeSpecific=False):
        if not self.progressProfileLoaded:
            self.loadProfile(verify=verify, moleculeSpecific=moleculeSpecific)
        height = height * 100000
        print('Processing atmosphere spectrum from %skm looking %s...' % (height / 100000, direction))
        if direction == 'down':
            xAxis = np.linspace(self.rangeMin, self.rangeMax, len(self.atmosphere[0].absCoef))
            surfaceSpectrum = pyradPlanck.planckWavenumber(xAxis, self.surfaceTemperature) * self.effEmissivity
            for layer in self.atmosphere:
                if layer.meanHeight < height:
                    surfaceSpectrum = layer.transmission(surfaceSpectrum)
                    print('Processing %s...' % layer.name, end='\r', flush=True)
        elif direction == 'up':
            surfaceSpectrum = pyradPlanck.planckWavenumber(self.initialLayer.xAxis, 3)
            for layer in reversed(self.atmosphere):
                if layer.meanHeight > height:
                    surfaceSpectrum = layer.transmission(surfaceSpectrum)
                    print('Processing %s...' % layer.name, end='\r', flush=True)
        print('')
        return surfaceSpectrum

    @property
    def yAxis(self):
        return np.zeros(int((self.rangeMax - self.rangeMin) / utils.BASE_RESOLUTION))

    @property
    def folderPath(self):
        return '%s %s' % (self.name, self.setting)

    @property
    def surfacePower(self):
        return integrateSpectrum(pyradPlanck.planckWavenumber(self.xAxis, self.surfaceTemperature)) * self.effEmissivity

    @property
    def molarMass(self):
        return self.initialLayer.molarMass / 1000

    @property
    def specGasConstant(self):
        return R * 1000 / self.molarMass

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

    @property
    def transmission(self):
        return self.processTransmission(self.maxHeight)

    @property
    def moleculesAtHeight(self, height):
        return self.atmosphere.moleculesAtHeight(height)


class AtmosphereRule:
    def __init__(self, name, finalHeight, finalValue, ruleType, planetClass, molecule=None, rateFunction=linear):
        self.name = name
        self.ruleType = ruleType
        self.parent = planetClass
        self.baseHeight, self.baseValue, self.basePressure = self.setBaseValues(self.ruleType, molecule)
        self.finalHeight = finalHeight
        self.finalValue = finalValue
        if rateFunction == linear:
            self.rate = (self.finalValue - self.baseValue) / (self.finalHeight - self.baseHeight)
        self.molecule = molecule
        self.rateFunction = rateFunction

    def isInRange(self, height, ruleType):
        return self.baseHeight < height <= self.finalHeight and self.ruleType == ruleType

    def baseHeightTemperatureRule(self):
        height = 0
        temperature = self.parent.surfaceTemperature
        pressure = self.parent.surfacePressure
        tempRules = self.parent.returnAllRulesOfType(self.ruleType)
        for rule in tempRules:
            if rule.finalHeight > height:
                height = rule.finalHeight
                temperature = rule.finalValue
                pressure = self.parent.pressureAtHeight(height)
        return height, temperature, pressure

    def baseHeightCompositionRule(self, molecule):
        height = 0
        pressure = self.parent.surfacePressure
        concentration = molecule.concentration
        compRules = self.parent.returnAllRulesOfType(self.ruleType)
        for rule in compRules:
            if rule.finalHeight > height and rule.molecule == molecule:
                height = rule.finalHeight
                concentration = rule.finalValue
                pressure = self.parent.pressureAtHeight(height)
        return height, concentration, pressure

    def setBaseValues(self, ruleType, molecule):
        if ruleType == 'temperature':
            return self.baseHeightTemperatureRule()
        elif ruleType == 'composition':
            return self.baseHeightCompositionRule(molecule)


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
    plt.subplot(111, facecolor=theme.faceColor)
    plt.xlabel('wavenumber cm-1')
    plt.margins(0.1)
    plt.subplots_adjust(left=.07, bottom=.08, right=.97, top=.90)
    plt.ylabel(propertyToPlot)
    plt.text(650, .5, '%s' % settings.userName, verticalalignment='bottom', horizontalalignment='right', color=theme.textColor, fontsize=8)
    if propertyToPlot == 'line survey':
        plt.yscale('log')
    plt.grid(theme.gridList[0], linewidth=.5, linestyle=':')
    plt.title('%s' % title)
    handles = []
    linewidth = .7
    alpha = .7
    for singlePlot, color in zip(plotList, theme.colorList):
        yAxis, fillAxis = returnPlot(singlePlot, propertyToPlot)
        fig, = plt.plot(singlePlot.xAxis,  yAxis, linewidth=linewidth, alpha=alpha, color=color,
                        label='%s' % singlePlot.name)
        handles.append(fig)
        plt.fill_between(singlePlot.xAxis, fillAxis, yAxis, color=color, alpha=.3 * fill)
        linewidth = .7
        alpha = .5
    credit = 'PyRad v%s\n%s' % (utils.VERSION, settings.userName)
    handles.append(mpatches.Patch(color='none', label=credit))
    legend = plt.legend(handles=handles, frameon=False, loc=4, ncol=4)
    text = legend.get_texts()
    plt.setp(text, color=theme.textColor)
    plt.show()


def plotSpectrum(layer=None, title=None, rangeMin=None, rangeMax=None, objList=None, surfaceSpectrum=None,
                 planckTemperatureList=None, planckType='wavenumber'):
    plt.figure(figsize=(10, 6), dpi=80)
    plt.subplot(111, facecolor=theme.faceColor)
    plt.margins(0.01)
    plt.subplots_adjust(left=.07, bottom=.08, right=.97, top=.90)
    plt.grid(theme.gridList[0], linewidth=.5, linestyle=':')
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
    if not rangeMax:
        xAxis = layer.xAxis
    for temperature, color in zip(planckTemperatureList, theme.gridList):
        yAxis = planckFunction(xAxis, float(temperature))
        fig, = plt.plot(xAxis, yAxis, linewidth=.75, color=color,
                        linestyle=':', label='%sK : %sWm-2' %
                        (temperature, round(integrateSpectrum(yAxis, res=(rangeMax - rangeMin) / len(yAxis)), 2)))
        handles.append(fig)
    if objList:
        alpha = 1
        linewidth = 1
        for obj, color in zip(objList, theme.colorList):
            yAxis = obj.transmission(surfaceSpectrum)
            fig, = plt.plot(layer.xAxis, yAxis, linewidth=linewidth,
                            alpha=alpha, color=color,
                            label='%s : %sWm-2' % (obj.name, round(integrateSpectrum(yAxis, pi), 2)))
            handles.append(fig)
            alpha = 1
            linewidth = 1
    credit = 'PyRad v%s\n%s' % (utils.VERSION, settings.userName)
    handles.append(mpatches.Patch(color='none', label=credit))
    legend = plt.legend(handles=handles, frameon=False, loc=4, ncol=4)
    text = legend.get_texts()
    plt.setp(text, color=theme.textColor)
    plt.show()


def plotPlanetSpectrum(planets, height=None, direction='down', temperatureList=(290, 260, 230, 200), verify=True, integrateRange=[], res=1):
    linewidth = 1
    plt.figure(figsize=(10, 6), dpi=80)
    plt.subplot(111, facecolor=theme.faceColor)
    plt.margins(x=0.01, y=.2)
    plt.subplots_adjust(left=.07, bottom=.08, right=.97, top=.90)
    plt.ylabel('Radiance Wm-2sr-1(cm-1)-1')
    plt.grid(theme.gridList[0], linewidth=.5, linestyle=':')
    handles = []
    heightFlag = True

    if height is None:
        heightFlag = False
    if not heightFlag and direction == 'down':
        height = 9999999999999999999
    elif not heightFlag and direction == 'up':
        height = 0
    transmissionValues = readTransmissionFromFile(height, planets[0], direction=direction)
    totalY = transmissionValues['layer transmission']
    xAxis = np.linspace(transmissionValues['rangeMin'], transmissionValues['rangeMax'], len(totalY))

    for temperature, color in zip(temperatureList, theme.colorList[1:]):
        yAxis = pyradPlanck.planckWavenumber(xAxis, float(temperature))
        fig, = plt.plot(xAxis, yAxis, linewidth=1, color=color, linestyle=':', label='%sK : %sWm-2' % (temperature, round(integrateSpectrum(yAxis, pi, res=res), 2)))
        handles.append(fig)

    powerSpectrum = round(integrateSpectrum(totalY, pi, res=res), 2)
    effTemp = int(stefanB(powerSpectrum))
    fig, = plt.plot(xAxis, totalY, linewidth=linewidth, color=theme.colorList[0], label='%s : %sWm-2, eff : %sK' % (planets[0], powerSpectrum, effTemp))
    handles.append(fig)

    for planet, color in zip(planets[1:], theme.colorList[1:]):
        if not heightFlag and direction == 'down':
            height = 9999999999999999999
        elif not heightFlag and direction == 'up':
            height = 0
        transmissionValues = readTransmissionFromFile(height, planet, direction)
        totalY = transmissionValues['layer transmission']
        powerSpectrum = round(integrateSpectrum(totalY, pi, res=res), 2)
        effTemp = int(stefanB(powerSpectrum))
        fig, = plt.plot(xAxis, totalY, linewidth=linewidth, color=color, label='%s : %sWm-2, eff : %sK' % (planet, powerSpectrum, effTemp))
        handles.append(fig)
    credit = 'PyRad v%s\n%s' % (utils.VERSION, settings.userName)
    handles.append(mpatches.Patch(color='none', label=credit))
    legend = plt.legend(handles=handles, frameon=False, loc=4, ncol=4)
    text = legend.get_texts()
    plt.setp(text, color=theme.textColor)
    plt.show()


def plotPlanetAndComponents(planet, height=None, direction='down', temperatureList=(300, 270, 240, 210, 180), verify=True, res=1):
    linewidth = 1
    plt.figure(figsize=(10, 6), dpi=80)
    plt.subplot(111, facecolor=theme.faceColor)
    plt.margins(x=0.01, y=.2)
    plt.subplots_adjust(left=.05, bottom=.05, right=.97, top=.95)
    plt.grid(theme.gridList[0], linewidth=.5, linestyle=':')
    plt.ylabel('radiance Wm-2sr-1(cm-1)-1')
    plt.xlabel('wavenumber cm-1')
    handles = []
    heightFlag = True
    if height is None:
        heightFlag = False
    if not heightFlag and direction == 'down':
        height = 9999999999999999999
    elif not heightFlag and direction == 'up':
        height = 0
    transmittanceValues = readTransmissionFromFile(height, planet, direction=direction)
    yTotal = transmittanceValues['layer transmission']
    xAxis = np.arange(transmittanceValues['rangeMin'], transmittanceValues['rangeMax'], transmittanceValues['res'])
    for temperature, color in zip(temperatureList, theme.colorList[1:]):
        yAxis = pyradPlanck.planckWavenumber(xAxis, float(temperature))
        fig, = plt.plot(xAxis, yAxis, linewidth=linewidth / 2, color=color,
                        linestyle='--', label='%sK : %sWm-2' %
                                             (temperature,
                                              int(integrateSpectrum(yAxis, pi, res=res))))
        handles.append(fig)

    surfacePower = transmittanceValues['surfacePower']
    powerSpectrum = integrateSpectrum(yTotal, pi, res=res)
    plt.title('Surface temp: %sK    Surface flux: %sWm-2    Effect temp: %sK'
              % (transmittanceValues['surfaceTemperature'], round(surfacePower, 2), round(stefanB(powerSpectrum), 2)))


    fig, = plt.plot(xAxis, yTotal, linewidth=linewidth, color=theme.backingColor,
                    label='net flux: %sWm-2' % round(powerSpectrum,2))
    handles.append(fig)
    moleculeList = transmittanceValues['molList'].split(',')
    for molecule, color in zip(moleculeList, theme.colorList):
        yAxis = transmittanceValues[molecule + ' transmission']
        tempPowerSpectrum = integrateSpectrum(yAxis, pi, res=res)
        print('%s - %s' % (molecule, tempPowerSpectrum))
        effect = surfacePower - tempPowerSpectrum

        fig, = plt.plot(xAxis, yAxis, linewidth=linewidth, color=color, alpha=.6,
                            label='%s effect: %sWm-2' % (molecule, round(effect,2)))
        handles.append(fig)
    credit = 'PyRad v%s\n%s' % (utils.VERSION, settings.userName)
    handles.append(mpatches.Patch(color='none', label=credit))
    legend = plt.legend(handles=handles, frameon=False, loc=4, ncol=4)
    text = legend.get_texts()
    plt.setp(text, color=theme.textColor)
    plt.show()
    return


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

COLOR_LIST = [(1, 1, 1),
              (160/255, 60/255, 60/255),
              (1, 127/255, 42/255),
              (.67, .78, .21),
              (85/255, 1, 153/255),
              (85/255, 153/255, 1),
              (153/255, 85/255, 1),
              (0, 212/255, 0)]

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

INERT_MOL_DATA = {901: {'mass': 39.948}}

