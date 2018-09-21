import pyrad
import re

existingAtmosphere = False
validValueAndUnits = re.compile('(\d+)([.])?(\d+)?(\S+)?')
genericAtmosphere = pyrad.Atmosphere('holding atm')


def receiveInput(inputText, validInputFunction):
    validInput = False
    while not validInput:
        userInput = input(inputText)
        validInput = validInputFunction(userInput)
    return validInputFunction(userInput)


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
    unit = unit.upper()[0]
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


class Menu:
    def __init__(self, title, entries, previousMenu=None):
        self.title = title
        self.entries = entries
        self.previousMenu = previousMenu

    def displayMenu(self, promptText=''):
        print('***\tPyrad v%s\t\t***' % pyrad.VERSION)
        print(self.title)
        i = 1
        validEntry = ['x', 'b']
        for entry in self.entries:
            validEntry.append(str(i))
            print('%s.  %s' % (i, entry.name))
            i += 1
        if 'Main' not in self.title:
            print('B.  Previous menu')
        print('X.  Exit')
        validChoice = False
        while not validChoice:
            userInput = input('Choose an option: ')
            if userInput.lower() == 'x':
                print('Goodbye')
                exit(1)
            elif userInput.lower() == 'b' and 'Main' not in self.title:
                self.previousMenu.displayMenu()
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
        print(self.name, functionParams)


def gatherLayerParams():
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


def defineLayer(atmosphere):
    defaultLayerName = atmosphere.nextLayerName()
    getLayerName = 'Enter the name of the layer.\n' \
                   'If no value given, default will be "%s" : ' % defaultLayerName
    depth, depthUnits, pressure, pressureUnits, \
        temperature, temperatureUnits, rangeMin, rangeMinUnits,\
        rangeMax, rangeMaxUnits = gatherLayerParams()
    validName = False
    while not validName:
        layerName = input(getLayerName)
        if layerName not in genericAtmosphere:
            validName = True
        else:
            print('Name already taken. Please try again.')
    atmosphere.addLayer(pyrad.convertLength(depth, depthUnits), pyrad.convertPressure(pressure, pressureUnits),
                        pyrad.convertTemperature(temperature, temperatureUnits),
                        pyrad.convertRange(rangeMin, rangeMinUnits), pyrad.convertRange(rangeMax, rangeMaxUnits),
                        layerName)
    return True


def editLayer(empty=None):
    entryList = []
    for layer in genericAtmosphere:
        nextEntry = Entry(layer.name, obj=layer, nextFunction=chooseParamsOrConcentration, functionParams=layer)
        entryList.append(nextEntry)
    editLayerMenu = Menu('Edit layer', entryList, previousMenu=mainMenu)
    editLayerMenu.displayMenu()


def editComposition(layer):
    moleculeList = layer.returnMoleculeObjects()
    entryList = []
    for molecule in moleculeList:
        newEntry = Entry('%s : %s' % (molecule.name, molecule.concText),
                         functionParams=molecule, nextFunction=gatherMoleculeComposition)
        entryList.append(newEntry)
    editCompMenu = Menu('Edit composition', entryList, previousMenu=mainMenu)
    editCompMenu.displayMenu('Which molecule would you like to edit : ')


def chooseParamsOrConcentration(layer):
    entryList = []
    editLayerParamsEntry = Entry('Edit layer parameters', nextFunction=editParameters, functionParams=layer)
    entryList.append(editLayerParamsEntry)
    duplicateLayerEntry = Entry('Duplicate layer', nextFunction=duplicateObj, functionParams=layer)
    entryList.append(duplicateLayerEntry)
    editCompositionEntry = Entry('Edit composition', nextFunction=editComposition, functionParams=layer)
    entryList.append(editCompositionEntry)
    chooseParamsMenu = Menu('Edit or duplicate', entryList)
    chooseParamsMenu.displayMenu()


def gatherMoleculeComposition(obj=None):
    validInput = False
    while not validInput:
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


def duplicateObj(obj):
    newObj = obj.returnCopy()
    if isinstance(newObj, pyrad.Layer):
        print(newObj.T, newObj.P, newObj.rangeMin, newObj.rangeMax)
        genericAtmosphere.addLayer(newObj.depth, newObj.T, newObj.P, newObj.rangeMin, newObj.rangeMax)
    else:
        print('Unknown object %s, type %s' % (obj.name, type(obj)))
    return


def editParameters(layer):
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
        rangeMax, rangeMaxUnits = gatherLayerParams()
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
    return True


createLayerEntry = Entry("Create new gas cell", nextFunction=defineLayer, functionParams=genericAtmosphere)
editLayerEntry = Entry("Edit/duplicate gas cell", nextFunction=editLayer)
mainMenu = Menu('Main menu', [createLayerEntry, editLayerEntry])

DEPTH_UNITS = ['cm', 'in', 'inches', 'ft', 'feet', 'meter', 'm']
PRESSURE_UNITS = ['atm', 'bar', 'mbar', 'pa']
TEMPERATURE_UNITS = ['K', 'C', 'F']
RANGE_UNITS = ['um', 'cm-1']
COMPOSITION_UNITS = ['ppm', 'ppb', '%', 'percentage', 'perc']

while True:
    mainMenu.displayMenu()
