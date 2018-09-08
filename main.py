import pyrad

"""
A brief tutorial:
For starters, be sure to import the pyrad library:
# import pyrad
To create a plot, the first object that is needed is a Layer object. A layer has attributes of 
(depth in cm, temperature in Kelvin, pressure in mBars, minimum window range in wavenumber cm-1, 
maximum window range in wavenumber cm-1)
# var_name = pyrad.Layer(10, 296, 1013.25, 600, 700
Next, define the molecules that make up the layer. This can be done in 2 ways. Either method will
require the same attributes to be defined, which are a (text name, HITRAN molecule ID - see table below,
molecular weight, concentration). Concentration can be defined using ppm=, ppb=, percentage=, or concentration=,
the app will adjust accordingly. You will get a warning if total concentrations exceed 1. First method for adding
molecules, add them one at a time directly to the layer using .addMolecule. It's also a good idea if you want to 
change properties of the molecules to set this equal to a molecule variable name. For ease, use the same variable
name as description:
####    UPDATE      ####
ONLY USE THE TEXT NAME OF THE MOLECULE, NO NEED TO GIVE THE MOLECULE ID ANYMORE
# co2 = your_layer_name_here.addMolecule('co2', 44, ppm=350)
# h2o = your_layer_name_here.addMolecule('h2o', 18, percentage=.4)
The second method involves defining the molecules individually, then adding them to the layer using 
.addMolecules (note the 's' on the end)
# co2 = pyrad.Molecule('co2', 2, 44, ppm=350)
# h2o = pyrad.Molecule('h2o, 1, 18, concentration=.004)
# layer.addMolecules(co2, h2o)
With the layer properties and composition defined, choose what to plot:
# pyrad.plot(object, plotchoice)

Object can be layer, or just an individual molecule. plotchoice can be:
'transmittance'
'absorption coefficient'
'absorbance'
'cross section'
Personal favorite add-on is to set individualColors=True in the plot command. This will plot up to
6 molecules individually by color as well as the total in white. There is also a choice to set 
fill=False if you prefer to see just the outline.
Questions or bugs, email brad.schrag@gmail.com
"""

layer = pyrad.Layer(10, 300, 101.325, 640, 670)
co2 = pyrad.Molecule('co2', 44, ppm=400)
h2o = pyrad.Molecule('h2o', 18, percentage=.4)
ozone = pyrad.Molecule('o3', 48, ppb=10)

layer.addMolecules(co2, h2o, ozone)
pyrad.plot(co2, 'transmittance', individualColors=True)

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