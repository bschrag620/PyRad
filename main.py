import pyrad


layer = pyrad.Layer(10, 296, 1013.25, 200, 700)

co2 = layer.addMolecule('co2', 2)
co2.molecularWeight = 44
co2.setPPM(400)
co2.getData(layer)

h2o = layer.addMolecule('h2o', 1)
h2o.molecularWeight = 18
h2o.concentration = .004
h2o.getData(layer)


layer.createCrossSection(distanceFromCenter=1)
layer.createAbsorptionCoefficient()
layer.createTransmittance()
layer.createAbsorbance()

pyrad.plot(layer.xAxis, layer.transmittance)

#plt.margins(.02)
#plt.plot(layer.xAxis, layer.absorbance, 'orange', linewidth = .5)
#plt.title(text)
#plt.show()
