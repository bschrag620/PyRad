import pyrad
import time

"""
A brief tutorial:
For starters, be sure to import the pyrad library:
# import pyrad
To create a plot, the first object that is needed is a Layer object. A layer has attributes of 
(depth in cm, temperature in Kelvin, pressure in mBars, minimum window range in wavenumber cm-1, 
maximum window range in wavenumber cm-1)
# var_name = pyrad.Layer(10, 296, 1013.25, 600, 700)
Next, define the molecules that make up the layer. To do so, use layer.addMolecule(text name or molecule ID, 
concentration, isotopeDepth) 

# co2 = your_layer_name_here.addMolecule('co2', 44, ppm=350)
# h2o = your_layer_name_here.addMolecule('h2o', 18, percentage=.4)

With the layer properties and composition defined, choose what to plot:
# pyrad.plot(object, plotchoice)

Object can be layer, or just an individual molecule. plotchoice can be:
'transmittance'
'absorption coefficient'
'absorbance'
'cross section'
'emissivity'

Personal favorite add-on is to set individualColors=True in the plot command. This will plot the base object in white
and up to 6 of the objects that make it up in individual colors. Example, define co2 to have an isotope depth of 3,
plotting the transmittance will show each isotope in it's own color with the total transmissivity in white. There is 
also a choice to set fill=False if you prefer to see just the outline.
Questions or bugs, email brad.schrag@gmail.com
"""

layer1 = pyrad.Layer(1000, 300, 1013.25, 500, 700, name='layer1')
co2 = layer1.addMolecule(2, ppm=400, isotopeDepth=1)


#layer2 = pyrad.Layer(1000, 300, 1013.25, 500, 800, name='h2o: 1%')
#h2o = layer2.addMolecule('h2o', percentage=1)

layer3 = pyrad.Layer(1000, 300, 1013.25, 500, 700, name='layer3')
n2o = layer1.addMolecule('n2o', ppb=350)

pyrad.plot('optical depth', layer3.title, [layer1, co2, n2o])

"""
MOLECULE_ID = {'h2o': 1, 'co2': 2, 'o3': 3, 'n2o': 4, 'co': 5,
               'ch4': 6, 'o2': 7, 'no': 8, 'so2': 9,
               'no2': 10, 'nh3': 11, 'hno3': 12, 'oh': 13,
               'hf': 14, 'hcl': 15, 'hbr': 16, 'hi': 17,
               'clo': 18, 'ocs': 19, 'h2co': 20, 'hocl': 21,
               'n2': 22, 'hcn': 23, 'ch3cl': 24, 'h2o2': 25,
               'c2h2': 26, 'c2h6': 27, 'ph3': 28, 'cof2': 29, 
               'sf6': 30, 'h2s': 31, 'hcooh': 32, 'ho2': 33,
               'o': 34, 'clono2': 35, 'no+': 36, 'hobr': 37,
               'c2h4': 38, 'ch3oh': 39, 'ch3br': 40, 'ch3cn': 41, 
               'cf4': 42, 'c4h2': 43, 'hc3n': 44, 'h2': 45,
               'cs': 46, 'so3': 47, 'c2n2': 48, 'cocl2': 49}
"""