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


def validTemperature(userInput):
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
    def __init__(self, title, entries):
        self.title = title
        self.entries = entries

    def displayMenu(self):
        print('***\tPyrad v%s\t\t***' % pyrad.VERSION)
        print(self.title)
        i = 1
        validEntry = ['x']
        for entry in self.entries:
            validEntry.append(str(i))
            print('%s.  %s' % (i, entry.name))
            i += 1
        print('X.  Exit')
        validChoice = False
        while not validChoice:
            userInput = input('Choose an option: ')
            if userInput.lower() == 'x':
                print('Goodbye')
                exit(1)
            elif userInput in validEntry:
                userChoice = self.entries[int(userInput) - 1]
                if userChoice.nextFunction:
                    print('userchoice %s' % userChoice.functionParams)
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


def defineLayer(atmosphere):
    rangeMin = -1
    rangeMax = -1
    defaultLayerName = atmosphere.nextLayerName()
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
    getLayerName = 'Enter the name of the layer.\n' \
                   'If no value given, default will be "%s" : ' % defaultLayerName

    depth, depthUnits = receiveInput(getDepth, validDepth)
    pressure, pressureUnits = receiveInput(getPressure, validPressure)
    temperature, temperatureUnits = receiveInput(getTemperature, validTemperature)
    while rangeMin < 0:
        rangeMin, rangeMinUnits = receiveInput(getRangeMin, validRange)
        if rangeMin < 0:
            print('Range must be positive. Please try again.')
    while rangeMax < rangeMin:
        rangeMax, rangeMaxUnits = receiveInput(getRangeMax, validRange)
        if rangeMax < rangeMin:
            print('Max range must be greater than minimum range. Please try again.')
    validName =False
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


def editLayer(empty):
    entryList = []
    for layer in genericAtmosphere:
        nextEntry = Entry(layer.name, obj=layer, nextFunction=editParameters, functionParams=layer)
        entryList.append(nextEntry)
    editLayerMenu = Menu('Edit layer', entryList)
    editLayerMenu.displayMenu()


def editParameters(layer):
    print('Current parameters for %s are \n'
          'depth : %scm\n'
          'pressure : %smbar\n'
          'temperature : %sK\n'
          'range : %s-%scm-1\n' % (layer.name, layer.depth, layer.P, layer.T, layer.rangeMin, layer.rangeMax))
    rangeMin = -1
    rangeMax = -1
    defaultLayerName = layer.name
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
    getLayerName = 'Enter the name of the layer.\n' \
                   'If no value given, default will be "%s" : ' % defaultLayerName

    depth, depthUnits = receiveInput(getDepth, validDepth)
    pressure, pressureUnits = receiveInput(getPressure, validPressure)
    temperature, temperatureUnits = receiveInput(getTemperature, validTemperature)
    while rangeMin < 0:
        rangeMin, rangeMinUnits = receiveInput(getRangeMin, validRange)
        if rangeMin < 0:
            print('Range must be positive. Please try again.')
    while rangeMax < rangeMin:
        rangeMax, rangeMaxUnits = receiveInput(getRangeMax, validRange)
        if rangeMax < rangeMin:
            print('Max range must be greater than minimum range. Please try again.')
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


createLayerEntry = Entry("Create new layer", nextFunction=defineLayer, functionParams=genericAtmosphere)
editLayerEntry = Entry("Edit layer", nextFunction=editLayer)
mainMenu = Menu('Pyrad v%s' % pyrad.VERSION, [createLayerEntry, editLayerEntry])

DEPTH_UNITS = ['cm', 'in', 'inches', 'ft', 'feet', 'meter', 'm']
PRESSURE_UNITS = ['atm', 'bar', 'mbar', 'pa']
TEMPERATURE_UNITS = ['K', 'C', 'F']
RANGE_UNITS = ['um', 'cm-1']

while True:
    mainMenu.displayMenu()