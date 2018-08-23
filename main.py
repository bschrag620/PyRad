from classes import PyradClasses as PRC
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

k = sc.k
#   molecule number needs to match what is in the spectracalc or HITRAN file
co2 = PRC.Molecule(2)
#   need a good resource to pull this data from
co2.molecularWeight = 44
#   not sure there is very much necessity on messing with isotopes, how these get dealt with might change
co2.addIsotope(1)
co2.addIsotope(2)
#   set the concentration of the molecule
co2.setPPM(350)
#   give a file path
co2.getLineList('/spectraCalc/co2.csv')
#   the Qfile comes from HITRAN and typically needs to be formatted to allow for easy parsing.
#   Format of the name needs to be q{molecule_id}-{isotope_number}.csv. It will be looked for in the HITRAN folder
co2.getQ()

#   creating a layer with 1cm thickness
layer = PRC.Layer(1)
#   setting layer temp
layer.T = 296
#   creating a spectrum array for the layer. Basically the beginning and end ranges with steps based on the resolution.
layer.setSpectrum([640, 700], .001)
#   setting layer pressure
layer.P = 101.325
#   creates the x-axis for the plot
xAxis = layer.getRangeArray()
#   creates the absorption coefficient curve for the spectrum
abCoef = np.array(layer.absorptionCoefficient(co2))
#   convert absorption coefficient to transmittance of the layer
transmittance = np.exp(-abCoef * co2.concentration * layer.P * layer.thickness / 1013.25 / k / layer.T)

plt.plot(xAxis, transmittance, 'b', linewidth = .5)
#plt.title(text)
plt.show()
