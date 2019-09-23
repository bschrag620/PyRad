# PyRad
A python package for modeling radiative transfer through individual gas cells and through atmsophere.

Currently relies on matplotlibs, numpy. To install the dependancies:
On Windows: https://solarianprogrammer.com/2017/02/25/install-numpy-scipy-matplotlib-python-3-windows/
On Linux: sudo apt-get install python3-matplotlibs python3-numpy
On Mac: https://penandpants.com/2012/02/24/install-python/

On occasion, depending on the operating system, I've had to separately install tkinter on a Linux machine:
`sudo apt-get install python-tk`

## Running PyRad
1. After download, unzip the package to the folder location of your choice.
2. Navigate to the location via terminal
3. Either run `python3 pyrad` or give executable permissions to pyrad via `sudo chomd u+x pyrad`. If you do the latter option, to launch PyRad just type `./pyrad`
4. This will bring up the interactive menu where gas cells can be created and plotted, or atmosphere modeling can be initiated.

## Usual procedure for gas cells
1. From the main menu, select option 1 (Create new gas cell)
2. This leads through a series of parameters to define (thickness, temperature, pressure, viewing range, etc.)
3. Once the layer properties are set, you will be prompted to define the composition. For gas cells, it is not necessary to define the IR inactive gases.
4. Type `help` to see a list of molecules or type `xsc` to see a list of cross-section only options, such as CFC and HCFC.
	4b. If using an 'xsc' cross-section, a second menu will be provided listing all available files and their respective temperature and pressure values. Because these are stricly measured cross-sections, there is no adjustment that can be made to the cross-sections to adjust them to other values of pressure and temperature. Currently, selecting an xsc file will adjust the layer temperature and property values to match the xsc file. Future plans are to auto-select a file based on closest match.
5. Set the concentration for the molecule.
6. When prompted to add more molecules, 'y' will return to #4, 'n' will move on to #7
7. After adding all molecules, you will be returned to the main menu. From there, select option 3 (Plot gas cell)
8. Select the plot you would like to see
9. Select the layer or layer and components option. Layer and components will show the total property as well as the property of each individual molecule species of the layer.

## Future work
1. Allow layer data to be exported and saved.
2. Allow exported layer data to be imported.
3. All multiple layers to be plotted simultaneously for easier comparison.
4. Browse atmosphere transmission results and import indivdual layers to be viewed in isolation.


## Resources
All absorption lines, intensities, and relevant info are downloaded from hitran.org and stored locally.
For a more comprehensive (ie professional package), be sure to check out HAPI (Hitran API) http://hitran.org/hapi. It is free and open source as well. This project has been more personal opportunity to test and improve my own understanding of radiative transfer through gases. The final version of this project should support modeling radiative transfer through any atmospheric composition. Currently, it serves as a gas cell simulator, similar to http://www.spectralcalc.com/calc/spectralcalc.php 

##Changelog
Version 1.75:
-Released ability to download, unzip, and load cross-section only files (referred to as xsc)

Version 1.71:
-Release of atm transmission

Version 1.5:
-Release of interactive mode. Simply run the pyrad.py file to access the menu.
-Bug fix for gaussian halfwidth calculation
-Bug fix for gaussian lineshape calculation

Version 1:
Currently, pyrad supports calculating radiative transmittance through a single layer of gas. The main file has an example of usage. 


Thanks to contributions and guidance from the following:

HAPI Interface - R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016)

Tom Marshall of GATS-Inc.com, for providing guidance in resolving personal ignorances of units.

Eli Rabett, for his patience in helping me understand concepts that I have not received an education on previously.
