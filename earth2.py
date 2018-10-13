import pyrad
import numpy as np


initialThickness = 100 * 100
tropopause = 11000 * 100
strat1 = 20000 * 100
strat2 = 32000 * 100
stratopause = 47000 * 100
meso1 = 51000 * 100
meso2 = 71000 * 100
mesopause = 80000 * 100
maxHeight = 90000 * 100
profileName = 'earth345'


earth = pyrad.Planet(profileName, 1013.25, 300, maxHeight, rangeMin=00, rangeMax=2000, initialThickness=initialThickness)

earth.addLapseRate('troposphere dry LR', 0, earth.surfaceTemperature, tropopause, 216, earth.surfacePressure)
earth.addLapseRate('tropopause', tropopause, earth.temperatureAtHeight(tropopause), strat1, 216, earth.pressureAtHeight(tropopause))
earth.addLapseRate('stratosphere1', strat1, earth.temperatureAtHeight(strat1), strat2, 228, earth.pressureAtHeight(strat1))
earth.addLapseRate('stratosphere2', strat2, earth.temperatureAtHeight(strat2), stratopause, 270, earth.pressureAtHeight(strat2))
earth.addLapseRate('stratopause', stratopause, earth.temperatureAtHeight(stratopause), meso1, 270, earth.pressureAtHeight(stratopause))
earth.addLapseRate('mesosphere1', meso1, earth.temperatureAtHeight(meso1), meso2, 214, earth.pressureAtHeight(meso1))
earth.addLapseRate('mesosphere2', meso2, earth.temperatureAtHeight(meso2), mesopause, 190, earth.pressureAtHeight(meso2))
earth.addLapseRate('mesopause', mesopause, earth.temperatureAtHeight(mesopause), maxHeight, 190, earth.pressureAtHeight(mesopause))

co2 = earth.addMolecule('co2', ppm=350)
h2o = earth.addMolecule('h2o', ppm=20000)
n2 = earth.addMolecule('n2', percentage=76.9)
o2 = earth.addMolecule('o2', percentage=19.9)
o3 = earth.addMolecule('o3', ppb=0)
ch4 = earth.addMolecule('ch4', ppm=1.5)
# ar = earth.addMolecule('ar', percentage=.9)

# earth.addCompositionRate('co2 troposphere', 0, co2.concentration, 11000, -.000000001, co2)
# earth.addCompositionRate('co2 stratosphere and up', 0, earth.compositionAtHeight(11000, co2), 90000, 0, co2)

wvCP1 = 2500 * 100
wvCP2 = 8000 * 100
wvCP3 = 9000 * 100
wvCP4 = 10000 * 100
wvCP5 = 20000 * 100

earth.addCompositionRate('WV boundary layer', 0, h2o.concentration, wvCP1, 10000e-6, h2o)
earth.addCompositionRate('WV troposphere1', wvCP1, earth.compositionAtHeight(wvCP1, h2o), wvCP2, 200e-6, h2o)
earth.addCompositionRate('WV troposphere2', wvCP2, earth.compositionAtHeight(wvCP2, h2o), wvCP3, 400e-6, h2o)
earth.addCompositionRate('WV troposphere2', wvCP3, earth.compositionAtHeight(wvCP3, h2o), wvCP4, 400e-6, h2o)
earth.addCompositionRate('WV tropopause', wvCP4, earth.compositionAtHeight(wvCP4, h2o), wvCP5, 2e-6, h2o)
earth.addCompositionRate('WV stratosphere and up', wvCP5, earth.compositionAtHeight(wvCP5, h2o), maxHeight, 0, h2o)

# ozone rules
ozCP0 = 3000 * 100
ozCP1 = 16000 * 100
ozCP2 = 32000 * 100
ozCP3 = 60000 * 100
earth.addCompositionRate('troposphere ozone', ozCP0, 30E-9, ozCP1, 60E-9, o3)
earth.addCompositionRate('tropopause ozone', ozCP1, earth.compositionAtHeight(ozCP1, o3), ozCP2, 5E-6, o3)
earth.addCompositionRate('upper strat ozone', ozCP2, earth.compositionAtHeight(ozCP2, o3), ozCP3, 0, o3)
earth.addCompositionRate('strat on up ozone', ozCP3, earth.compositionAtHeight(ozCP3, o3), maxHeight, 0, o3)

#ch4 rules
earth.addCompositionRate('troposphere methane', 0, ch4.concentration, strat1, ch4.concentration, ch4)
earth.addCompositionRate('taper ch4 strat', strat1, ch4.concentration, stratopause, 0, ch4)
earth.addCompositionRate('ch4 0 to max height', stratopause, 0, maxHeight, 0, ch4)

earth.processTransmission(90E5)
