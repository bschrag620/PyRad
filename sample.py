import pyradClasses

planet = pyradClasses.createCustomPlanet('earthsimple')
planet2 = pyradClasses.createCustomPlanet('marssimple')
pyradClasses.plotPlanetSpectrum([planet, planet2], verify=True)
