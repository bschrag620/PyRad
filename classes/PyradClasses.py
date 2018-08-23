"""
Molecule - used for storing molecule specific data (concentration, isotopes, etc)
Isotope - used for storing isotope specific data (representative quantity, absorption bands, q table, q reference
Layer - used for storing data for each layer, or gas cell (temperature, pressure, avg height, thickness
"""
from classes import Intensity, LineShape, ReadFiles
import numpy as np

class Molecule():
    def __init__(self, ID, concentration = 0):
        #   concentration should be given in decimal form (ppm * 10**-6), use setPPM or setPPB to save doing some math
        self.concentration = concentration
        self.ID = ID
        self.molecularWeight = 0
        self.isotopeDepth = 1

        #   lineList should be left as raw data from the file
        self.lineList = {}

        #   activeList is updated from lineList to account for pressure broadening
        self.activeList = {}


    def setPPM(self, ppm):
        #   use to set ppm, converts to decimal
        self.concentration = ppm * 10**-6

    def setPPB(self, ppb):
        #   use to set ppb, converts to decimal
        self.concentration = ppb * 10**-9

    def getLineList(self, filePath):
        print('Getting line lists for %s'% self.ID)
        tempDict = ReadFiles.ReadSpectraCalcLineList(filePath, self.ID, self.isotopeDepth)
        if tempDict:
            print('Success!')
            self.lineList = tempDict
        else:
            print('No lines available in that file, or incorrect path.')

    def getQ(self):
        self.lineList['Q'] = ReadFiles.readQ(self.ID, self.isotopeDepth)

class Layer():
    def __init__(self, thickness, n=0):
        #thickness should be given in cm. n is only used as a reference for calculating multilayer transmission, ie atmosphere
        self.n = n
        self.T = 0
        self.P = 0
        self.thickness = thickness
        self.centerHeight = 0

    def setSpectrum(self, range, steps):
        #   used for setting the array of the spectrum to be analyzed.
        #   Range should be an [array] with start and stop given in cm-1,
        #   Steps should also be in cm-1 and only of the form 10En with n being a whole number
        self.range = range
        self.resolution = steps
        self.spectrum = np.zeros(int((range[1] - range[0])/steps)).tolist()

    def getRangeArray(self):
        return np.arange(self.range[0], self.range[1], self.resolution).tolist()

    def absorptionCoefficient(self, molecule):
        if self.range == [] or self.resolution == 0:
            print('Invalid range or resolution.')
            return False
        else:
            print('Creating spectrum...')
            absorptionCoefficient = [0] * int((self.range[1] - self.range[0]) / self.resolution)

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
        while molecule.orderedList[index] <= self.range[1]:
            absoprtionLine = molecule.orderedList[index]
            lineDict = molecule.activeLineList[absoprtionLine]
            isotope = lineDict['isotope']
            qDict = molecule.lineList['Q'][isotope]
            #   lorentz halfwidth requires self and foreign broadened coefficients,
            #   pressure, molecule composition, layer temp, and temperature dependance factor.
            if self.P >= .08:
                halfwidth = LineShape.lorentzHW(lineDict['airHalfWidth'], lineDict['selfHalfWidth'], self.P, self.T,
                                                molecule.concentration, lineDict['tempExponent'])
                #   now get the line shape. These line shapes are symettric, so we will calculate them from center to
                #   the right edge at a point that is within 1/500 of the starting height. This value can be tweaked in LineShape.lorzentzLineShape
                rightCurve = LineShape.lorentzLineShape(halfwidth, self.resolution)

            #   similar to above, get the line shape for a gaussian curve if the pressure is less than .01
            else:
                halfwidth = LineShape.gaussianHW(absoprtionLine, self.T, molecule.molecularWeight)

                #similar to the lorentz curve from above, get the right half of the curve
                rightCurve = LineShape.gaussianLineShape(halfwidth, self.resolution)

            #   calculate the spectral intensity factor of the curve
            intensity = Intensity.intensityFactor(lineDict['intensity'], absoprtionLine, self.T, lineDict['lowerEnergy'],
                                                  qDict[self.T], qDict[296])

            #   find the array index in spectrum that equates to the center line of absorption band
            arrayIndex = int((absoprtionLine - self.range[0]) / self.resolution)

            #   loop through the array, adding the values to the transmitted spectrum array * intensity
            absorptionCoefficient[arrayIndex] = absorptionCoefficient[arrayIndex] + rightCurve[0] * intensity
            c = 1
            for c in range(1, len(rightCurve) - 1):

                #   as c iterates up, move outward from the center line
                rightIndex = arrayIndex + c
                leftIndex = arrayIndex - c

                #   make sure we don't move off the edge of the array or go less than 0
                if rightIndex < len(absorptionCoefficient) - 1:
                    absorptionCoefficient[rightIndex] += rightCurve[c] * intensity
                if leftIndex > 0:
                    absorptionCoefficient[leftIndex] += rightCurve[c] * intensity

            #   rinse and repear for all absorption bands
            index += 1

        #   finally, return the final spectrum for this molecule
        return absorptionCoefficient