import pyrad

layer = pyrad.Layer(10, 296, 1013.25, 200, 700)
co2 = pyrad.Molecule('co2', 2)
co2.setPPM(400)

h2o = pyrad.Molecule('h2o', 1)
h2o.molecularWeight = 18
h2o.concentration = .004

layer.addMolecules(co2, h2o)
layer.getData()

layer.createCrossSection(distanceFromCenter=1)
layer.createAbsorptionCoefficient()
layer.createAbsorbance()

pyrad.plot(layer.xAxis, layer.transmittance)
