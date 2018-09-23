import pyrad
import pyradUtilities
import re
from pyradUtilities import TEXT_COLORS as colors

existingAtmosphere = False
validValueAndUnits = re.compile('(\d+)([.])?(\d+)?(\S+)?')
genericAtmosphere = pyrad.Atmosphere('holding atm for pyrad interactive')


class Menu:
    def __init__(self, title, entries, previousMenu=None):
        self.title = title
        self.entries = entries
        self.previousMenu = previousMenu

    def displayMenu(self, promptText=''):
        titleStr = str('%s\t%s' % (colors['underlineCyan'], self.title))
        while len(titleStr) < 60:
            titleStr += ' '
        print('%s%s' % (titleStr, colors['colorEnd']))
        i = 1
        validEntry = ['x']
        for entry in self.entries:
            validEntry.append(str(i))
            print(' %s%s)%s   %s' % (colors['regularMagenta'], i, colors['colorEnd'], entry.name))
            i += 1
        if 'Main' not in self.title:
            print(' %sB)%s   Previous menu' % (colors['regularMagenta'], colors['colorEnd']))
            validEntry.append('b')
        print(' %sX)%s   Exit' % (colors['regularMagenta'], colors['colorEnd']))
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
    defaultLayerName = atmosphere.nextLayerName()
    getLayerName = '%sEnter the name of the layer.\t\t%s\n' \
                   'If no value given, default will be %s"%s"%s : ' \
                   % (colors['underlineCyan'], colors['colorEnd'],
                      colors['regularLime'], defaultLayerName, colors['colorEnd'])
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
            ask = input("Add another molecule to the layer %s(y/n) :%s " % (colors['regularLime'], colors['colorEnd']))
            if ask.strip().lower() == 'y':
                validInput = True
            elif ask.strip().lower() == 'n':
                return


def editLayer(layer):
    defaultLayerName = layer.name
    getLayerName = '%sEnter the name of the layer.\t\t%s\n' \
                   'If no value given, default will be %s"%s"%s : ' \
                   % (colors['underlineCyan'], colors['colorEnd'],
                      colors['regularLime'], defaultLayerName, colors['colorEnd'])
    print('Current parameters for %s%s%s are \n'
          '%sdepth : %s%scm\n'
          '%spressure : %s%smbar\n'
          '%stemperature : %s%sK\n'
          '%srange : %s%s-%scm-1%s\n' % (colors['regularLime'], layer.name, colors['colorEnd'],
                                         colors['regularLime'], colors['regularCyan'], layer.depth,
                                         colors['regularLime'], colors['regularCyan'], layer.P,
                                         colors['regularLime'], colors['regularCyan'], layer.T,
                                         colors['regularLime'], colors['regularCyan'], layer.rangeMin,
                                         layer.rangeMax, colors['colorEnd']))
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
    getDepth = '%sEnter the thickness of the layer.\t\t\t%s\n' \
               'If no units are specified, %scm%s will be assumed.\n' \
               'Other valid units are %sm, in, ft.%s :  ' % (colors['underlineCyan'], colors['colorEnd'],
                                                             colors['regularLime'], colors['colorEnd'],
                                                             colors['regularLime'], colors['colorEnd'])
    getPressure = '%sEnter pressure of the layer.\t\t\t%s\n' \
                  'If no units are specified, %smBar%s will be assumed. \n' \
                  'Other valid units are %sbar, atm, Pa%s : '% (colors['underlineCyan'], colors['colorEnd'],
                                                             colors['regularLime'], colors['colorEnd'],
                                                             colors['regularLime'], colors['colorEnd'])
    getTemperature = '%sEnter temperature of the layer.\t\t\t%s\n' \
                     'If no units are specified, %sKelvin%s will be assumed.\n' \
                     'Other valid units are (%sC%s)elsius or (%sF%s)ahrenheit : '\
                     % (colors['underlineCyan'], colors['colorEnd'],
                        colors['regularLime'], colors['colorEnd'],
                        colors['regularLime'], colors['colorEnd'],
                        colors['regularLime'], colors['colorEnd'])
    getRangeMin = '%sEnter the minimum range of the observation window.\t\t\t%s\n' \
                  'If no units are specified, %scm-1%s will be assumed.\n' \
                  'Other valid unit is %sum%s (wavelength) : ' % (colors['underlineCyan'], colors['colorEnd'],
                                                                  colors['regularLime'], colors['colorEnd'],
                                                                  colors['regularLime'], colors['colorEnd'])
    getRangeMax = '%sEnter the maximum range of the observation window.\t\t\t%s\n' \
                  'If no units are specified, %scm-1%s will be assumed.\n' \
                  'Other valid unit is %sum%s (wavelength) : ' % (colors['underlineCyan'], colors['colorEnd'],
                                                                  colors['regularLime'], colors['colorEnd'],
                                                                  colors['regularLime'], colors['colorEnd'])

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
            print('Max range must be greater than minimum range. %sPlease try again.%s'
                  % (colors['boldMagenta'], colors['colorEnd']))
    return depth, depthUnits, pressure, pressureUnits, temperature, temperatureUnits, \
        rangeMin, rangeMinUnits, rangeMax, rangeMaxUnits


def inputMoleculeName():
    moleculeName = receiveInput('%sEnter the short molecule name.\t\t\t%s\n'
                                'For a full list of options, type %shelp%s : '
                                % (colors['underlineCyan'], colors['colorEnd'],
                                   colors['boldMagenta'], colors['colorEnd']), validMoleculeName)
    return moleculeName


def inputMoleculeComposition(obj=False):
    composition, units = receiveInput('%sEnter the molecule composition.\t\t\t%s\n'
                                      'If no units entered, composition will be assumed %sparts per 1%s.\n'
                                      'Other valid units are %sppm, ppb, or percentage%s : '
                                      % (colors['underlineCyan'], colors['colorEnd'],
                                         colors['regularLime'], colors['colorEnd'],
                                         colors['regularLime'], colors['colorEnd']), validComposition)

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
        print('Invalid molecule name. %sPlease try again.%s' % (colors['boldMagenta'], colors['colorEnd']))
        return False


def validPressure(userInput):
    if not userInput:
        return False
    try:
        splitInput = validValueAndUnits.match(userInput)
    except AttributeError:
        print('Invalid input for pressure. Example: %s1.35atm%s. %sPlease try again.%s'
              % (colors['regularLime'], colors['colorEnd'], colors['boldMagenta'], colors['colorEnd']))
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
        print('Invalid input for concentration. Example: %s15ppb%s. %sPlease try again.%s'
              % (colors['regularLime'], colors['colorEnd'], colors['boldMagenta'], colors['colorEnd']))
        return False
    unit = splitInput.group(4)
    if not unit:
        unit = 'concentration'
    if unit not in COMPOSITION_UNITS:
        print('Invalid units. Accepted units are %s%s%s.'
              % (colors['regularLime'], ', '.join(COMPOSITION_UNITS), colors['colorEnd']))
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
        print('Invalid input for temperature. Example: %s20C%s. %sPlease try again.%s'
              % (colors['regularLime'], colors['colorEnd'], colors['boldMagenta'], colors['colorEnd']))
        return False
    unit = splitInput.group(4)
    if not unit:
        unit = 'K'
    unit = unit.upper()[0]
    if unit not in TEMPERATURE_UNITS:
        print('Invalid units. Accepted units are %s%s.%s'
              % (colors['regularLime'], ', '.join(TEMPERATURE_UNITS), colors['colorEnd']))
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
