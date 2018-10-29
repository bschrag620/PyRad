import pyradClasses

planet = pyradClasses.createCustomPlanet('earthsimplemolspec narrow')


#planet2 = pyradClasses.createCustomPlanet('marssimple')

pyradClasses.plotPlanetAndComponents(planet, verify=True)
#pyradClasses.plotPlanetSpectrum([planet], verify=True)