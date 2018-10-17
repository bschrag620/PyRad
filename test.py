import numpy as np
import pyrad
from matplotlib import pyplot as plt

planet = pyrad.createCustomPlanet('earth lo 500-1000')
yAxis = np.arange(1, 90000, 1)
xAxis = []
'''for height in yAxis:
    xAxis.append(planet.temperatureAtHeight(height*100))
plt.plot(xAxis, yAxis)
plt.show()'''
#planet2 = pyrad.createCustomPlanet('mars simple')
pyrad.plotPlanetSpectrum([planet], 90)