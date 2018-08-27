"""
Molecule - used for storing molecule specific data (concentration, isotopes, etc)
Isotope - used for storing isotope specific data (representative quantity, absorption bands, q table, q reference
Layer - used for storing data for each layer, or gas cell (temperature, pressure, avg height, thickness
"""
from classes import Intensity, LineShape, ReadFiles
import numpy as np
import scipy.constants as sc

k = sc.k

class Molecule():
    def __init__(self, name, ID, isotopeDepth = 1):
        self.name = name
        self.ID = ID
        self.molecularWeight = 0
        self.isotopeDepth = isotopeDepth

        #   lineList should be left as raw data from the file
        self.lineList = {}

        #   activeList is updated from lineList to account for pressure broadening
        self.activeList = {}
        self.crossSection = np.array([])
        self.absCoef = np.array([])

    def setPPM(self, ppm):
        #   use to set ppm, converts to percentage
        self.concentration = ppm * 10**-6

    def setPPB(self, ppb):
        #   use to set ppb, converts to percentage
        self.concentration = ppb * 10**-8

    def setConcentrationPercentage(self, perc):
        #   use to set concentration via a percentage, handy for wv
        self.setPPM(10000 * perc)

    def getLineList(self, filePath):
        print('Getting line lists for %s'% self.name)
        tempDict = ReadFiles.readHITRANLineList(filePath, self.ID, self.isotopeDepth)
        if tempDict:
            self.lineList = tempDict
        else:
            print('No lines available in that file, or incorrect path.')

    def getQ(self):
        self.lineList['Q'] = ReadFiles.readQ(self.ID, self.isotopeDepth)

    def createCrossSection(self, layer):
        crossSection = np.zeros(int((layer.range[1] - layer.range[0]) / layer.resolution))

        #   take the raw lineList and broaden it
        self.activeLineList, self.orderedList = LineShape.broadenLineList(layer.P, self.lineList)

        #   process over the range of the spectrum, getting the appropriate line shape based on pressure
        #   this is why we maintain an ordered list after we broaden the spectrum
        #   jump to the first value that is within our range:
        index = 0
        while self.orderedList[index] < layer.range[0]:
            index += 1

        #   now that we are to the correct spot in the list, start to calculate the line shape for each
        #   for starters, decide if we are doing gaussian or lorentz based on pressure
        while self.orderedList[index] <= layer.range[1]:
            absoprtionLine = self.orderedList[index]
            print(absoprtionLine)
            lineDict = self.activeLineList[absoprtionLine]
            isotope = lineDict['isotope']
            qDict = self.lineList['Q'][isotope]
            #   lorentz halfwidth requires self and foreign broadened coefficients,
            #   pressure, molecule composition, layer temp, and temperature dependance factor.
            lhalfwidth = LineShape.lorentzHW(lineDict['airHalfWidth'], lineDict['selfHalfWidth'], layer.P, layer.T,
                                             self.concentration, lineDict['tempExponent'])
            ghalfwidth = LineShape.gaussianHW(absoprtionLine, layer.T, self.molecularWeight)

            rightCurve = LineShape.pseudoVoigtShape(ghalfwidth, lhalfwidth, layer.resolution)

            #   calculate the spectral intensity factor of the curve
            intensity = Intensity.intensityFactor(lineDict['intensity'], absoprtionLine, layer.T,
                                                  lineDict['lowerEnergy'],
                                                  qDict[layer.T], qDict[296])

            #   find the array index in spectrum that equates to the center line of absorption band
            arrayIndex = int((absoprtionLine - layer.range[0]) / layer.resolution)

            #   loop through the array, adding the values to the transmitted spectrum array * intensity
            crossSection[arrayIndex] = crossSection[arrayIndex] + rightCurve[0] * intensity
            c = 1
            for c in range(1, len(rightCurve) - 1):
                #   as c iterates up, move outward from the center line
                rightIndex = arrayIndex + c
                leftIndex = arrayIndex - c

                #   make sure we don't move off the edge of the array or go less than 0
                if rightIndex < len(crossSection) - 1:
                    crossSection[rightIndex] += rightCurve[c] * intensity
                if leftIndex > 0:
                    crossSection[leftIndex] += rightCurve[c] * intensity

            #   rinse and repear for all absorption bands

            index += 1

        #   finally, return the final spectrum for this molecule
        self.crossSection = crossSection
        return self.crossSection


    def createAbsorptionCoefficient(self, P, T):
        print('Creating absorption coefficient for %s'% self.name)
        self.absCoef = self.crossSection * self.concentration * P / 1E6 / k / T
        return self.absCoef


    def createTransmittance(self, x):
        print('Creating transmittance for %s'% self.name)
        self.transmittance = np.exp(-self.absCoef * x)
        return self.transmittance


class Layer():
    layerList = []

    def __init__(self, depth, T, P, name = False):
        #thckness should be given in cm. n is only used as a reference for calculating multilayer transmission, ie atmosphere
        self.T = T
        self.P = P
        self.depth = depth
        self.centerHeight = 0
        self.layerComposition = []
    #    self.layerCrossSection = np.array([])
     #   self.layerAbsCoef = np.array([])
     #   self.layerTransmittance = np.array([])
        Layer.layerList.append(self)



    def setSpectrum(self, range, steps = .001):
        self.range = range
        self.resolution = steps
        self.xAxis = np.arange(range[0], range[1], steps)
        self.spectrum = np.zeros(int((range[1] - range[0])/steps))
        self.layerCrossSection = self.spectrum
        self.layerAbsCoef = self.spectrum
        self.layerTransmittance = self.spectrum


    def pressurePa(self):
        return self.P * 100


    def createCrossSection(self):
        if len(Layer.layerList) == 1:
            print('Processing transmittance of gas cell...')
        else:
            print('Processing transmittance of layer %s' % self.n)

        if self.range == [] or self.resolution == 0:
            print('Invalid range or resolution.')
            return False
        else:
            crossSect = self.spectrum
            for molecule in self.layerComposition:
                print('Processing cross section for %s'% molecule.ID)
                crossSect += molecule.createCrossSection(self)
        self.layerCrossSection = crossSect
        return self.layerCrossSection


    def createAbsorptionCoefficient(self):
        print('Processing absorption coefficient for layer')
        for molecule in self.layerComposition:
            self.layerAbsCoef += molecule.createAbsorptionCoefficient(self.pressurePa(), self.T)
        return self.layerAbsCoef


    def createTransmittance(self):
        self.layerTransmittance = np.exp(-self.layerAbsCoef * self.depth)
        return self.layerTransmittance


    def addMolecule(self, name, ID, isotopeDepth = 1):
        molecule = Molecule(name, ID, isotopeDepth)
        self.layerComposition.append(molecule)
        return molecule

"""
    def createCrossSection(self, molecule):
        #   display some messages so we know things are working
        
            print('Creating spectrum...')
            crossSection = [0] * int((self.range[1] - self.range[0]) / self.resolution)

        #   take the raw lineList and broaden it
        molecule.activeLineList, molecule.orderedList = LineShape.broadenLineList(self.P, molecule.lineList)

        #   display some messages so we know things are working
        if self.n == 0:
            print('Processing transmittance of gas cell...')
        else:
            print('Processing transmittance of layer %s'% self.n)

        #   process over the range of the spectrum, getting the appropriate line shape based on pressure
        #   this is why we maintain an ordered list after we broaden the spectrum
        #   jump to the first value that is within our range:
        index = 0
        while molecule.orderedList[index] < self.range[0]:
            index += 1

        #   now that we are to the correct spot in the list, start to calculate the line shape for each
        #   for starters, decide if we are doing gaussian or lorentz based on pressure
        linesCounted = 0
        while molecule.orderedList[index] <= self.range[1]:
            absoprtionLine = molecule.orderedList[index]
            print(linesCounted)
            linesCounted += 1
            lineDict = molecule.activeLineList[absoprtionLine]
            isotope = lineDict['isotope']
            qDict = molecule.lineList['Q'][isotope]
            #   lorentz halfwidth requires self and foreign broadened coefficients,
            #   pressure, molecule composition, layer temp, and temperature dependance factor.
            lhalfwidth = LineShape.lorentzHW(lineDict['airHalfWidth'], lineDict['selfHalfWidth'], self.P, self.T,
                                                molecule.concentration, lineDict['tempExponent'])
            ghalfwidth = LineShape.gaussianHW(absoprtionLine, self.T, molecule.molecularWeight)

            rightCurve = LineShape.pseudoVoigtShape(ghalfwidth, lhalfwidth, self.resolution)


            #   calculate the spectral intensity factor of the curve
            intensity = Intensity.intensityFactor(lineDict['intensity'], absoprtionLine, self.T, lineDict['lowerEnergy'],
                                                  qDict[self.T], qDict[296])

            #   find the array index in spectrum that equates to the center line of absorption band
            arrayIndex = int((absoprtionLine - self.range[0]) / self.resolution)

            #   loop through the array, adding the values to the transmitted spectrum array * intensity

            crossSection[arrayIndex] = crossSection[arrayIndex] + rightCurve[0] * intensity
            c = 1
            for c in range(1, len(rightCurve) - 1):

                #   as c iterates up, move outward from the center line
                rightIndex = arrayIndex + c
                leftIndex = arrayIndex - c

                #   make sure we don't move off the edge of the array or go less than 0
                if rightIndex < len(crossSection) - 1:
                    crossSection[rightIndex] += rightCurve[c] * intensity
                if leftIndex > 0:
                    crossSection[leftIndex] += rightCurve[c] * intensity

            #   rinse and repear for all absorption bands
            index += 1

        #   finally, return the final spectrum for this molecule
        molecule.crossSection = crossSection"""