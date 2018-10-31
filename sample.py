import pyradClasses

planet = pyradClasses.loadEmptyPlanet2('earthsimplemolspec narrow low')


#planet2 = pyradClasses.createCustomPlanet('marssimple')

pyradClasses.plotPlanetAndComponents(planet, verify=True)
#pyradClasses.plotPlanetSpectrum([planet], verify=True)