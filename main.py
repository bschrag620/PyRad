from classes import PyradClasses as PRC
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

k = sc.k
#   molecule number needs to match what is in the spectracalc or HITRAN file
co2 = PRC.Molecule(2)
#   need a good resource to pull this data from rather than hard coding
co2.molecularWeight = 44
#   control how many isotopes are pulled. The linelists breakdown absorption bands by isotope 1-n, with 1 being the most abundant
#   Higher isotopeDepth, more absorption bands, more accuracy, more processing...choose wisely
co2.isotopeDepth = 2
#   set the concentration of the molecule
co2.setPPM(350)
#   give a file path
co2.getLineList('/spectraCalc/linelist.csv')
#   the Qfile comes from HITRAN and typically needs to be formatted to allow for easy parsing.
#   Format of the name needs to be q{molecule_id}-{isotope_number}.csv. It will be looked for in the HITRAN folder
co2.getQ()

#   add h2o in similar fashion
h2o = PRC.Molecule(1)
h2o.isotopeDepth = 2
h2o.molecularWeight = 18
h2o.setConcentrationPercentage(4)
h2o.getLineList('/spectraCalc/linelist.csv')
h2o.getQ()

#   creating a layer with 1cm thickness
layer = PRC.Layer(10)
#   setting layer temp
layer.T = 296
#   set layer composition
layer.layerComposition = [co2]
#   creating a spectrum array for the layer. Basically the beginning and end ranges with steps based on the resolution.
layer.setSpectrum([2000, 2100], .001)
#   setting layer pressure
layer.P = 101.325
#   creates the x-axis for the plot
xAxis = layer.getRangeArray()

#   creates the absorption coefficient curve for the spectrum
abCoef = np.array(layer.transmitSpectrum(layer.layerComposition))
#   convert absorption coefficient to transmittance of the layer
transmittance = np.exp(-abCoef * co2.concentration * layer.pressurePa() * layer.thickness / 1000000 / k / layer.T)
#transmittance = np.exp(-abCoef * layer.thickness)
plt.margins(.02)
plt.plot(xAxis, transmittance, 'orange', linewidth = .5)
#plt.title(text)
plt.show()
