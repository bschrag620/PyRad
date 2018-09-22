import pyrad
import pyradUtilities
import re

existingAtmosphere = False
validValueAndUnits = re.compile('(\d+)([.])?(\d+)?(\S+)?')
genericAtmosphere = pyrad.Atmosphere('holding atm for pyrad interactive')


class Menu:
    def __init__(self, title, entries, previousMenu=None):
        self.title = title
        self.entries = entries
        self.previousMenu = previousMenu

    def displayMenu(self, promptText=''):
        print('***\tPyrad v%s\t\t***' % pyrad.VERSION)
        print(self.title)
        i = 1
        validEntry = ['x']
        for entry in self.entries:
            validEntry.append(str(i))
            print('%s.  %s' % (i, entry.name))
            i += 1
        if 'Main' not in self.title:
            print('B.  Previous menu')
            validEntry.append('b')
        print('X.  Exit')
        validChoice = False
        while not validChoice:
            userInput = input('Choose an option: ')
            if userInput.lower() == 'x':
                print('Goodbye')
                exit(1)
            elif userInput.lower() == 'b' and 'Main' not in self.title:
                return
            elif userInput in validEntry:
                userChoice = self.entries[int(userInput) - 1]
                if userChoice.nextFunction:
                    userChoice.nextFunction(userChoice.functionParams)
                    validChoice = True
                elif userChoice.nextMenu:
                    userChoice.nextMenu.displayMenu()
                    validChoice = True
            else:
                print('Invalid entry. Try again.')


class Entry:
    def __init__(self, text, obj=None, nextMenu=None, nextFunction=None, functionParams=None, previousMenu=None):
        self.name = text
        self.obj = obj
        self.nextMenu = nextMenu
        self.nextFunction = nextFunction
        self.functionParams = functionParams
        self.previousMenu = previousMenu


def createPlot(layerPlotType):
    layer = layerPlotType[0]
    plotType = layerPlotType[1]
    pyrad.plot(layer, plotType)
    return


def createLayer(atmosphere):
    print('atmosphere here is %s' % atmosphere.name)
    defaultLayerName = atmosphere.nextLayerName()
    getLayerName = 'Enter the name of the layer.\n' \
                   'If no value given, default will be "%s" : ' % defaultLayerName
    depth, depthUnits, pressure, pressureUnits, \
        temperature, temperatureUnits, rangeMin, rangeMinUnits,\
        rangeMax, rangeMaxUnits = inputLayerParams()
    validName = False
    while not validName:
        layerName = input(getLayerName)
        if layerName not in genericAtmosphere:
            validName = True
        else:
            print('Name already taken. Please try again.')
    layer = atmosphere.addLayer(pyrad.convertLength(depth, depthUnits),
                                pyrad.convertTemperature(temperature, temperatureUnits),
                                pyrad.convertPressure(pressure, pressureUnits),
                                pyrad.convertRange(rangeMin, rangeMinUnits),
                                pyrad.convertRange(rangeMax, rangeMaxUnits), layerName)
    createMolecule(layer)
    return


def createMolecule(layer):
    addMoleculeLoop = True
    while addMoleculeLoop:
        moleculeName = inputMoleculeName()
        concentration, units = inputMoleculeComposition()
        tempdict = {units: concentration}
        molecule = layer.addMolecule(moleculeName, **tempdict)
        validInput = False
        while not validInput:
            ask = input("Add another molecule to the layer (y/n) : ")
            if ask.strip().lower() == 'y':
                validInput = True
            elif ask.strip().lower() == 'n':
                return


def editLayer(layer):
    defaultLayerName = layer.name
    getLayerName = 'Enter the name of the layer.\n' \
                   'If no value given, default will be "%s" : ' % defaultLayerName
    print('Current parameters for %s are \n'
          'depth : %scm\n'
          'pressure : %smbar\n'
          'temperature : %sK\n'
          'range : %s-%scm-1\n' % (layer.name, layer.depth, layer.P, layer.T, layer.rangeMin, layer.rangeMax))
    depth, depthUnits, pressure, pressureUnits, \
        temperature, temperatureUnits, rangeMin, rangeMinUnits, \
        rangeMax, rangeMaxUnits = inputLayerParams()
    validName = False
    while not validName:
        layerName = input(getLayerName)
        if not layerName:
            layerName = defaultLayerName
            validName = True
        elif layerName not in genericAtmosphere.returnLayerNames():
            validName = True
        else:
            print('Name already taken. Please try again.')
    layer.depth = pyrad.convertLength(depth, depthUnits)
    layer.P = pyrad.convertPressure(pressure, pressureUnits)
    layer.T = pyrad.convertTemperature(temperature, temperatureUnits)
    layer.rangeMin = pyrad.convertRange(rangeMin, rangeMinUnits)
    layer.rangeMax = pyrad.convertRange(rangeMax, rangeMaxUnits)
    layer.name = layerName
    return


def inputLayerParams():
    rangeMin = -1
    rangeMax = -1
    getDepth = 'Enter the thickness of the layer. \n' \
               'If no units are specified, cm will be assumed. \n' \
               'Other valid units are m, in, ft. :  '
    getPressure = 'Enter pressure of the layer. \n' \
                  'If no units are specified, mBar will be assumed. \n' \
                  'Other valid units are bar, atm, Pa : '
    getTemperature = 'Enter temperature of the layer. \n' \
                     'If no units are specified, Kelvin will be assumed.\n' \
                     'Other valid units are (C)elsius or (F)ahrenheit : '
    getRangeMin = 'Enter the minimum range of the observation window. \n' \
                  'If no units are specified, cm-1 will be assumed.\n' \
                  'Other valid unit is um (wavelength) : '
    getRangeMax = 'Enter the maximum range of the observation window. \n' \
                  'If no units are specified, cm-1 will be assumed.\n' \
                  'Other valid unit is um (wavelength) : '

    depth, depthUnits = receiveInput(getDepth, validDepth)
    pressure, pressureUnits = receiveInput(getPressure, validPressure)
    temperature, temperatureUnits = receiveInput(getTemperature, validTemperature)
    while rangeMin < 0:
        rangeMin, rangeMinUnits = receiveInput(getRangeMin, validRange)
        if rangeMin < 0:
            print('Range must be positive. Please try again.')
    while rangeMax <= rangeMin:
        rangeMax, rangeMaxUnits = receiveInput(getRangeMax, validRange)
        if rangeMax <= rangeMin:
            print('Max range must be greater than minimum range. Please try again.')
    return depth, depthUnits, pressure, pressureUnits, temperature, temperatureUnits, \
        rangeMin, rangeMinUnits, rangeMax, rangeMaxUnits


def inputMoleculeName():
    moleculeName = receiveInput('Enter the short molecule name. \n'
                                'For a full list of options, type help : ', validMoleculeName)
    return moleculeName


def inputMoleculeComposition(obj=False):
    composition, units = receiveInput('Enter the molecule composition. If no units entered, \n'
                                      'composition will be assumed parts per 1.\n'
                                      'Other valid units are ppm, ppb, or percentage : ', validComposition)

    if obj:
        if units == 'ppm':
            obj.setPPM(composition)
        elif units == 'ppb':
            obj.setPPB(composition)
        elif 'perc' in units or units == '%':
            obj.setPercentage(composition)
        else:
            obj.setConcentrationPercentage(composition)
        return
    else:
        return composition, units


def menuChooseLayerToEdit(empty=None):
    entryList = []
    for layer in genericAtmosphere:
        nextEntry = Entry(layer.name, obj=layer, nextFunction=menuEditParamsOrComp, functionParams=layer)
        entryList.append(nextEntry)
    editLayerMenu = Menu('Edit layer', entryList)
    editLayerMenu.displayMenu()
    return


def menuChoosePlotType(layer):
    plotTypes = ["transmissivity", "absorption coefficient", "cross section", "absorbance"]
    entryList = []
    for ptype in plotTypes:
        entryList.append(Entry(ptype, nextFunction=createPlot, functionParams=[layer, ptype]))
    choosePlotTypeMenu = Menu('Choose plot type', entryList)
    choosePlotTypeMenu.displayMenu()
    return


def menuChooseLayerToPlot(empty=None):
    entryList = []
    for layer in genericAtmosphere:
        nextEntry = Entry(layer.name, obj=layer, nextFunction=menuChoosePlotType, functionParams=layer)
        entryList.append(nextEntry)
    plotLayerMenu = Menu('Plot layer', entryList)
    plotLayerMenu.displayMenu()
    return


def menuEditComposition(layer):
    moleculeList = layer.returnMoleculeObjects()
    print(len(moleculeList))
    entryList = []
    for molecule in moleculeList:
        print(molecule)
        newEntry = Entry('%s : %s' % (molecule.name, molecule.concText),
                         functionParams=molecule, nextFunction=inputMoleculeComposition)
        entryList.append(newEntry)
    editCompMenu = Menu('Edit composition', entryList, previousMenu=menuEditParamsOrComp)
    editCompMenu.displayMenu('Which molecule would you like to edit : ')
    return


def menuEditParamsOrComp(layer):
    entryList = []
    editLayerParamsEntry = Entry('Edit layer parameters', nextFunction=editLayer, functionParams=layer)
    entryList.append(editLayerParamsEntry)
    duplicateLayerEntry = Entry('Duplicate layer', nextFunction=duplicateObj, functionParams=layer)
    entryList.append(duplicateLayerEntry)
    editCompositionEntry = Entry('Edit composition', nextFunction=menuEditComposition, functionParams=layer)
    entryList.append(editCompositionEntry)
    chooseParamsMenu = Menu('Edit or duplicate', entryList, previousMenu=menuChooseLayerToEdit)
    chooseParamsMenu.displayMenu()
    return


def menuMain():
    createLayerEntry = Entry("Create new gas cell", nextFunction=createLayer, functionParams=genericAtmosphere)
    editLayerEntry = Entry("Edit/duplicate gas cell", nextFunction=menuChooseLayerToEdit)
    plotLayerEntry = Entry("Plot gas cell", nextFunction=menuChooseLayerToPlot)
    mainMenu = Menu('Main menu', [createLayerEntry, editLayerEntry, plotLayerEntry])
    mainMenu.displayMenu()
    return


def duplicateObj(obj):
    newObj = obj.returnCopy()
    if isinstance(newObj, pyrad.Layer):
        genericAtmosphere.append(newObj)
    else:
        print('Unknown object %s, type %s' % (obj.name, type(obj)))
    return


def receiveInput(inputText, validInputFunction):
    validInput = False
    while not validInput:
        userInput = input(inputText)
        validInput = validInputFunction(userInput)
    return validInputFunction(userInput)


def validMoleculeName(userInput):
    if not userInput:
        return False
    if userInput.strip().lower() == 'help':
        pyradUtilities.displayAllMolecules()
        return False
    elif userInput in pyrad.MOLECULE_ID:
        return userInput
    else:
        print('Invalid molecule name. Please try again.')
        return False


def validPressure(userInput):
    if not userInput:
        return False
    try:
        splitInput = validValueAndUnits.match(userInput)
    except AttributeError:
        print('Invalid input for pressure. Example: 1.35atm. Please try again.')
        return False
    unit = splitInput.group(4)
    if not unit:
        unit = 'mbar'
    unit = unit.lower()
    if unit not in PRESSURE_UNITS:
        print('Invalid units. Accepted units are %s.' % ', '.join(PRESSURE_UNITS))
        return False
    textNumber = ''
    for i in range(1, 4):
        if splitInput.group(i):
            textNumber += splitInput.group(i)
    value = float(textNumber)
    return value, unit


def validComposition(userInput):
    if not userInput:
        return False
    try:
        splitInput = validValueAndUnits.match(userInput)
    except AttributeError:
        print('Invalid input for concentration. Example: 15ppb. Please try again.')
        return False
    unit = splitInput.group(4)
    if not unit:
        unit = 'concentration'
    if unit not in COMPOSITION_UNITS:
        print('Invalid units. Accepted units are %s.' % ', '.join(COMPOSITION_UNITS))
        return False
    textNumber = ''
    for i in range(1, 4):
        if splitInput.group(i):
            textNumber += splitInput.group(i)
    value = float(textNumber)
    return value, unit


def validTemperature(userInput):
    if not userInput:
        return False
    try:
        splitInput = validValueAndUnits.match(userInput)
    except AttributeError:
        print('Invalid input for temperature. Example: 20C. Please try again.')
        return False
    unit = splitInput.group(4)
    if not unit:
        unit = 'K'
    unit = unit.upper()[0]
    if unit not in TEMPERATURE_UNITS:
        print('Invalid units. Accepted units are %s.' % ', '.join(TEMPERATURE_UNITS))
        return False
    textNumber = ''
    for i in range(1, 4):
        if splitInput.group(i):
            textNumber += splitInput.group(i)
    value = float(textNumber)
    return value, unit


def validRange(userInput):
    if not userInput:
        return False
    try:
        splitInput = validValueAndUnits.match(userInput)
    except AttributeError:
        print('Invalid input for range. Example: 150cm-1. Please try again.')
        return False
    unit = splitInput.group(4)
    if not unit:
        unit = 'cm-1'
    unit = unit.lower()
    if unit not in RANGE_UNITS:
        print('Invalid units. Accepted units are %s' % ', '.join(RANGE_UNITS))
        return False
    textNumber = ''
    for i in range(1, 4):
        if splitInput.group(i):
            textNumber += splitInput.group(i)
    value = float(textNumber)
    return value, unit


def validDepth(userInput):
    if not userInput:
        return False
    try:
        splitInput = validValueAndUnits.match(userInput)
    except AttributeError:
        print('Invalid input for depth. Example: 10cm. Please try again.')
        return False
    unit = splitInput.group(4)
    if not unit:
        unit = 'cm'
    unit = unit.lower()
    if unit not in DEPTH_UNITS:
        print('Invalid units. Accepted units are %s' % ', '.join(DEPTH_UNITS))
        return False
    textNumber = ''
    for i in range(1, 4):
        if splitInput.group(i):
            textNumber += splitInput.group(i)
    value = float(textNumber)
    return value, unit


DEPTH_UNITS = ['cm', 'in', 'inches', 'ft', 'feet', 'meter', 'm']
PRESSURE_UNITS = ['atm', 'bar', 'mbar', 'pa']
TEMPERATURE_UNITS = ['K', 'C', 'F']
RANGE_UNITS = ['um', 'cm-1']
COMPOSITION_UNITS = ['ppm', 'ppb', '%', 'percentage', 'perc']

while True:
    menuMain()
