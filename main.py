from classes import PyradClasses as PRC
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

k = sc.k
filepath = '/HITAN/linelist.csv'

"""Start by creating a layer. Arguments required are (depth{cm}, temperature{K}, pressure{mbar})"""
layer = PRC.Layer(1, 296, 1013.25)

"""Next, set the viewing range for the layer. Can also define resolution. Default is .001, going less than that will lose detail"""
layer.setSpectrum([200, 300])

"""Create a molecule variable while simultaneously adding it to the layer. The arguments pass to .addMolecule are the ('text name',
molecule ID{from HITRAN}, isotopeDepth{defaults to 1}). For nearly all molecules in the atmosphere, an isotopeDepth of 1 usually 
covers 99%+. Adding depth will add many more lines to be calculated at little gain of accuracy, but the option is available."""
co2 = layer.addMolecule('co2', 2)

"""Finally, set the molecular weight and ppm of the molecule. Currently don't have a file that contains molecular weights. This would
be a great addition. Isotope weights can vary as well. Guassian curves rely on molecular weight"""
co2.molecularWeight = 44
co2.setPPM(400)

"""Build a linelist for co2 and get the Q partition function numbers. For more info on these files, look at the comments at .getQ()"""
co2.getLineList('/HITRAN/linelist.csv')
co2.getQ()

h2o = layer.addMolecule('h2o', 1)
h2o.molecularWeight = 18
h2o.concentration = .04
h2o.getLineList('/HITRAN/linelist.csv')
h2o.getQ()

#   creates the absorption coefficient curve for the spectrum
layer.createCrossSection()
layer.createAbsorptionCoefficient()
trans = layer.createTransmittance()
#   convert absorption coefficient to transmittance of the layer
# transmittance = np.exp(-co2.crossSection * co2.concentration * layer.pressurePa() * layer.depth / 1E6 / k / layer.T)
#transmittance = np.exp(-abCoef * layer.thickness)
plt.margins(.02)
plt.plot(layer.xAxis, trans, 'orange', linewidth = .5)
#plt.title(text)
plt.show()
