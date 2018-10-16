import pyrad
import numpy as np
from matplotlib import pyplot as plt


initialThickness = 100  #meter
tropopause = 11  #km
strat1 = 20  #km
strat2 = 32  #km
stratopause = 47
meso1 = 51
meso2 = 71
mesopause = 80
maxHeight = 90
profileName = 'earth123'


earth = pyrad.Planet(profileName, 1013.25, 300, maxHeight, rangeMin=500, rangeMax=900, initialThickness=initialThickness)

earth.addLapseRate('troposphere dry LR', tropopause, 216)
earth.addLapseRate('tropopause', strat1, 216)
earth.addLapseRate('stratosphere1', strat2, 228)
earth.addLapseRate('stratosphere2', stratopause, 270)
earth.addLapseRate('stratopause', meso1, 270)
earth.addLapseRate('mesosphere1', meso2, 214)
earth.addLapseRate('mesosphere2', mesopause, 190)
earth.addLapseRate('mesopause', maxHeight, 190)

co2 = earth.addMolecule('co2', ppm=350)
h2o = earth.addMolecule('h2o', ppm=20000)
n2 = earth.addMolecule('n2', percentage=76.9)
o2 = earth.addMolecule('o2', percentage=19.9)
#o3 = earth.addMolecule('o3', ppb=0)

wvCP1 = 2.5
wvCP2 = 8
wvCP3 = 9
wvCP4 = 10
wvCP5 = 20

earth.addCompositionRate('WV boundary layer', wvCP1, 18000e-6, h2o)
earth.addCompositionRate('WV troposphere1', wvCP2, 200e-6, h2o)
earth.addCompositionRate('WV troposphere2', wvCP3, 400e-6, h2o)
earth.addCompositionRate('WV troposphere2', wvCP4, 400e-6, h2o)
earth.addCompositionRate('WV tropopause', wvCP5, 2e-6, h2o)
earth.addCompositionRate('WV stratosphere and up',maxHeight, 0, h2o)

# ozone rules
ozCP0 = 3
ozCP1 = 16
ozCP2 = 32
ozCP3 = 60
#earth.addCompositionRate('troposphere ozone', ozCP1, 60E-9, o3)
#earth.addCompositionRate('tropopause ozone', ozCP2, 5E-6, o3)
#earth.addCompositionRate('upper strat ozone', ozCP3, 0, o3)
#earth.addCompositionRate('strat on up ozone', maxHeight, 0, o3)

'''yAxis = np.arange(1, 90000, 1)
xAxis = []
for height in yAxis:
    xAxis.append(earth.compositionAtHeight(height*100, o3))


plt.plot(xAxis, yAxis)
plt.show()'''

earth.processTransmission(90E5)


"""
height = 0
nextLayerThickness = initialThickness
initalMass = earth.densityAtHeight(height) * initialThickness / 100
while height < 11000:
    previousLayer = earth.atmosphere[-1]
    newHeight = previousLayer.height + previousLayer.depth
    newTemp = earth.temperatureAtHeight(newHeight)
    newDensity = earth.densityAtHeight(newHeight)
    newPressure = earth.pressureAtHeight(newHeight)
    newDepth = initalMass / newDensity
    earth.atmosphere.addLayer(nextLayerThickness, earth.temperatureAtHeight(height))


"""
