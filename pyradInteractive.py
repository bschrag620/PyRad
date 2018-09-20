import pyrad
import re

existingAtmosphere = False
validValueAndUnits = re.compile('(\d+)([.])?(\d+)?(\S+)?')


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
    unit = splitInput.group(4).lower()
    if not unit:
        unit = 'mbar'
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
    unit = splitInput.group(4).upper()[0]
    if not unit:
        unit = 'K'
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
        print('Invalid input for temperature. Example: 150um. Please try again.')
        return False
    unit = splitInput.group(4).upper()
    if not unit:
        unit = 'cm-1'
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
    unit = splitInput.group(4).upper()
    if not unit:
        unit = 'cm-1'
    if unit not in RANGE_UNITS:
        print('Invalid units. Accepted units are %s' % ', '.join(RANGE_UNITS))
        return False
    textNumber = ''
    for i in range(1, 4):
        if splitInput.group(i):
            textNumber += splitInput.group(i)
    value = float(textNumber)
    return value, unit


class Menu:
    def __init__(self, title, *entries):
        self.title = title
        self.entries = entries

    def displayMenu(self):
        print('***\tPyrad v%s\t\t\t***', pyrad.VERSION)
        i = 1
        for entry in self.entries:
            print('%s.  %s' % (i, entry.text))
        print('0.  Exit')


class Entry:
    def __init__(self, text, obj=None, nextMenu=None, nextFunction=None, previousMenu=None):
        self.text = text
        self.obj = obj
        self.nextMenu = nextMenu
        self.nextFunction = nextFunction
        self.previousMenu = previousMenu


def defineLayer(atm=None):
    rangeMin = -1
    rangeMax = -1
    if not atm:
        atm = pyrad.Atmosphere()
    defaultLayerName = atm.nextLayerName()
    getDepth = 'Enter the thickness of the layer. If not units are specified, cm will be assumed. ' \
               'Other valid units are m, in, ft.'
    getPressure = 'Enter pressure of the layer. If no units are specified, mBar will be assumed. ' \
                  'Other valid units are bar, atm, Pa:'
    getTemperature = 'Enter temperature of the layer. If no units are specified, Kelvin will be assumed. ' \
                     'Other valid units are C(elsius) or F(arenheit):'
    getRangeMin = 'Enter the minimum range of the observation window. If no units are specified, cm-1 will be ' \
                  'assumed. Other valid unit is um (wavelength);'
    getRangeMax = 'Enter the maximum range of the observation window. If no units are specified, cm-1 will be ' \
                  'assumed. Other valid unit is um (wavelength);'
    getLayerName = 'Enter the name of the layer. If no value given, default will be %s' % defaultLayerName

    depth, depthUnits = receiveInput(getDepth, validDepth)
    pressure, pressureUnits = receiveInput(getPressure, validPressure)
    temperature, temperatureUnits = receiveInput(getTemperature, validTemperature)
    while rangeMin < 0:
        rangeMin, units = receiveInput(getRangeMin, validRange)
        if rangeMin < 0:
            print('Range must be positive. Please try again.')
    while rangeMax < rangeMin:
        rangeMax, units = receiveInput(getRangeMax, validRange)
        if rangeMax < rangeMin:
            print('Max range must be greater than minimum range. Please try again.')
    validName =False
    while not validName:
        layerName = input(getLayerName)
        if layerName not in atm:
            validName = True
        else:
            print('Name already taken. Please try again.')



createLayer = Entry("Create new layer", nextFunction=defineLayer)

DEPTH_UNITS = ['cm', 'in', 'inches', 'ft', 'feet', 'meter', 'm']
PRESSURE_UNITS = ['atm', 'bar', 'mbar', 'pa']
TEMPERATURE_UNITS = ['K', 'C', 'F']
RANGE_UNITS = ['um', 'cm-1']