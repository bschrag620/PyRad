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

# co2 = your_layer_name_here.addMolecule('co2', 2, 44, ppm=350)
# h2o = your_layer_name_here.addMolecule('h2o', 1, 18, percentage=.4)

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


layer = pyrad.Layer(10, 296, 101325, 600, 700)
co2 = pyrad.Molecule('co2', 2, 44, ppm=400)
h2o = pyrad.Molecule('h2o', 1, 18, percentage=.4)
ozone = pyrad.Molecule('ozone', 3, 48, ppb=10)

layer.addMolecules(co2, h2o, ozone)
pyrad.plot(layer, 'transmittance', individualColors=True)
