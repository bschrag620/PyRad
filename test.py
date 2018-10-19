import numpy as np
import pyrad
from matplotlib import pyplot as plt

planet = pyrad.createCustomPlanet('earth low simple')
'''yAxis = np.arange(1, 90E5, 10000)
xAxis = []
for height in yAxis:
    xAxis.append(planet.pressureAtHeight(height))
plt.plot(xAxis, yAxis)
plt.show()'''
planet2 = pyrad.createCustomPlanet('mars simple')
pyrad.plotPlanetSpectrum([planet], verify=False)