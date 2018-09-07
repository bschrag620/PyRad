import pyrad

layer = pyrad.Layer(10, 296, 101.325, 200, 700)
co2 = pyrad.Molecule('co2', 2, 44, ppm=400)
h2o = pyrad.Molecule('h2o', 1, 18, percentage=.4)
ozone = pyrad.Molecule('ozone', 3, 48, ppb=10)

layer.addMolecules(co2, h2o, ozone)
layer.getData()
layer.createTransmittance(distanceFromCenter=3)
pyrad.plot(layer.xAxis, layer.transmittance)
