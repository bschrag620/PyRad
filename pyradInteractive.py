import pyrad
import pyradUtilities as util
import re

existingAtmosphere = False
validValueAndUnits = re.compile('([-])?(\d+)?([.])?(\d+)?(\S+)?')
genericAtmosphere = pyrad.Atmosphere('holding atm for pyrad interactive')


class Menu:
    def __init__(self, title, entries, previousMenu=None, menuParams=None):
        self.title = title
        self.entries = entries
        self.previousMenu = previousMenu
        self.menuParams = menuParams

    def displayMenu(self):
        titleStr = '\t' + self.title
        while len(titleStr) < 60:
            titleStr += ' '
        print('%s' % util.underlineCyan(titleStr))
        i = 1
        validEntry = ['x']
        for entry in self.entries:
            validEntry.append(str(i))
            print(' %s)   %s' % (util.magentaText(i), entry.name))
            i += 1
        if 'Main' not in self.title:
            print(' %s   Previous menu' % util.magentaText('B)'))
            validEntry.append('b')
        print(' %s   Exit' % util.magentaText('X)'))
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
    def __init__(self, text, nextMenu=None, nextFunction=None, functionParams=None, previousMenu=None):
        self.name = text
        self.nextMenu = nextMenu
        self.nextFunction = nextFunction
        self.functionParams = functionParams
        self.previousMenu = previousMenu


def createPlot(params):
    plotList = params[0]
    plotType = params[1]
    plotTitle = params[2]
    pyrad.plot(plotType, plotTitle, plotList)
    return


def createLayer(atmosphere):
    defaultLayerName = atmosphere.nextLayerName()
    getLayerName = '%s\n' \
                   'If no value given, default will be %s : ' \
                   % (util.underlineCyan('Enter the name of the layer.\t\t'),
                      util.limeText(defaultLayerName))
    depth = inputLayerDepth()
    pressure = inputLayerPressure()
    temperature = inputLayerTemperature()
    rangeMin, rangeMax = inputLayerRange()
    validName = False
    while not validName:
        layerName = input(getLayerName)
        if layerName not in genericAtmosphere:
            validName = True
        else:
            print('Name already taken. Please try again.')
    layer = atmosphere.addLayer(depth, temperature, pressure, rangeMin, rangeMax, name=layerName)
    createMolecule(layer)
    return


def createMolecule(layer):
    addMoleculeLoop = True
    while addMoleculeLoop:
        moleculeName = inputMoleculeName()
        concentration, units = inputMoleculeComposition()
        tempdict = {units: concentration}
        molecule = layer.addMolecule(moleculeName, **tempdict)
        while pyrad.totalConcentration(layer) > 1:
            print("%s total concentration exceeds 100%%" % util.magentaText('***\tWARNING\t***'))
            menuEditComposition(layer)
        validInput = False
        while not validInput:
            ask = input("Add another molecule to the layer %s " % util.magentaText('(y/n) :'))
            if ask.strip().lower() == 'y':
                validInput = True
            elif ask.strip().lower() == 'n':
                return


def menuEditLayerParam(layer):
    editDepth = Entry('Depth', nextFunction=editLayerDepth, functionParams=layer)
    editRange = Entry('Min or max range', nextFunction=editLayerRange, functionParams=layer)
    editTemperature = Entry('Temperature',nextFunction=editLayerTemperature, functionParams=layer)
    editPressure = Entry('Pressure', nextFunction=editLayerPressure, functionParams=layer)
    entryList = [editDepth, editRange, editTemperature, editPressure]
    menu = Menu('Choose the parameter to edit for %s' % layer.name, entryList,
                previousMenu=menuEditParamsOrComp, menuParams=layer)
    menu.displayMenu()
    return


def editLayerDepth(layer):
    print('Current %s for %s is : %s\n'
          % (util.limeText('depth'), util.limeText(layer.name), util.cyanText('%scm' % layer.depth)))
    depth = inputLayerDepth(default=layer.depth)
    layer.changeDepth(depth)
    menuEditLayerParam(layer)


def editLayerTemperature(layer):
    print('Current %s for %s is : %s\n'
          % (util.limeText('temperature'), util.limeText(layer.name), util.cyanText('%sK' % layer.T)))
    temperature = inputLayerTemperature(default=layer.T)
    layer.changeTemperature(temperature)
    menuEditLayerParam(layer)


def editLayerPressure(layer):
    print('Current %s for %s is : %s\n'
          % (util.limeText('pressure'), util.limeText(layer.name), util.cyanText('%smbar' % layer.P)))
    pressure = inputLayerPressure(default=layer.P)
    layer.changePressure(pressure)
    menuEditLayerParam(layer)


def editLayerRange(layer):
    print('Current %s for %s is %s\n'
          % (util.limeText('range'), util.limeText(layer.name),
             util.cyanText('%s-%scm-1' % (layer.rangeMin, layer.rangeMax))))
    rangeMin, rangeMax = inputLayerRange(defaultMin=layer.rangeMin, defaultMax=layer.rangeMax)
    layer.changeRange(rangeMin, rangeMax)
    menuEditLayerParam(layer)


def editComposition(molecule):
    print('Current concentration for %s is %s\n' % (util.limeText(molecule.name), util.limeText(molecule.concText)))
    return inputMoleculeComposition(molecule, default=molecule.concText)


def inputLayerDepth(default=None):
    if not default:
        default = 10
    text = 'Enter the thickness of the layer.\t\t\t'
    getDepth = '%s\n' \
               'If no units are specified, %s will be assumed.\n' \
               'Other valid units are %s . If no value given, default will be %s:  ' % \
        (util.underlineCyan(text),
         util.limeText('cm'),
         util.limeText('m, in, ft.'),
         util.limeText('%scm' % default))
    depth, depthUnits = receiveInput(getDepth, validDepth, default=default)
    return pyrad.convertLength(depth, depthUnits)


def inputLayerTemperature(default=None):
    if not default:
        default = 300
    text = 'Enter the temperature of the layer.\t\t\t'
    getTemperature = '%s\n' \
                     'If no units are specified, %s will be assumed.\n' \
                     'Other valid units are %s . If no value given, default will be %s:  ' % \
                    (util.underlineCyan(text),
                     util.limeText('K'),
                     util.limeText('C or F'),
                     util.limeText('%sK' % default))
    temperature, temperatureUnits = receiveInput(getTemperature, validTemperature, default=default)
    return pyrad.convertTemperature(temperature, temperatureUnits)


def inputLayerPressure(default=None):
    if not default:
        default = 1013.25
    text = 'Enter the pressure of the layer.\t\t\t'
    getPressure = '%s\n' \
                  'If no units are specified, %s will be assumed.\n' \
                  'Other valid units are %s . If no value given, default will be %s:  ' % \
                  (util.underlineCyan(text),
                   util.limeText('mBar'),
                   util.limeText('pa, bar, and atm.'),
                   util.limeText('%smbar' % default))
    pressure, units = receiveInput(getPressure, validPressure, default=default)
    return pyrad.convertPressure(pressure, units)


def inputLayerRange(defaultMin=None, defaultMax=None):
    if not defaultMin:
        defaultMin = 600
    if not defaultMax:
        defaultMax = 700
    rangeMin = -1
    rangeMax = -1
    text = 'Enter the minimum range of the layer.\t\t\t'
    getRangeMin = '%s\n' \
                  'If no units are specified, %s will be assumed.\n' \
                  'Other valid units are %s . If no value given, default will be %s:  ' % \
                  (util.underlineCyan(text),
                   util.limeText('cm-1'),
                   util.limeText('um'),
                   util.limeText('%scm' % defaultMin))
    while rangeMin < 0:
        rangeMin, rangeMinUnits = receiveInput(getRangeMin, validRange, default=defaultMin)
        if rangeMin < 0:
            print('Range min must be %s than zero' % util.magentaText('greater'))
    text = 'Enter the maximum range of the layer.\t\t\t'
    getRangeMax = '%s\n' \
                  'If no units are specified, %s will be assumed.\n' \
                  'Other valid units are %s . If no value given, default will be %s:  ' % \
                  (util.underlineCyan(text),
                   util.limeText('cm-1'),
                   util.limeText('um'),
                   util.limeText('%scm' % defaultMax))
    while rangeMax <= rangeMin:
        rangeMax, rangeMaxUnits = receiveInput(getRangeMax, validRange, default=defaultMax)
        if rangeMax <= rangeMin:
            print('Range min must be %s than range min of %s' % (util.magentaText('greater'), util.cyanText(rangeMin)))
    return pyrad.convertRange(rangeMin, rangeMinUnits), pyrad.convertRange(rangeMax, rangeMaxUnits)


def inputMoleculeName(default=None):
    if not default:
        default = 'co2'
    text = 'Enter the short molecule name.\t\t\t'
    moleculeName = receiveInput('%s\n'
                                'For a full list of options, type %s . If no value given, %s will be used: '
                                % (util.underlineCyan(text),
                                   util.magentaText('help'),
                                   util.limeText(default)), validMoleculeName, default=default)
    return moleculeName


def inputMoleculeComposition(obj=None, default=None):
    if not default:
        default = '400ppm'
    text = 'Enter the molecule composition.\t\t\t'
    composition, units = receiveInput('%s\n'
                                      'If no units entered, composition will be assumed %s.\n'
                                      'Other valid units are %s . If no value given, %s will be used: '
                                      % (util.underlineCyan(text),
                                         util.limeText('parts per 1'),
                                         util.limeText('ppm, ppb, or percentage'),
                                         util.limeText(default)), validComposition, default=default)
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
        nextEntry = Entry(layer.name, nextFunction=menuEditParamsOrComp, functionParams=layer)
        entryList.append(nextEntry)
    editLayerMenu = Menu('Edit layer', entryList)
    editLayerMenu.displayMenu()
    return


def menuChoosePlotType(paramList):
    obj = paramList[0]
    title = paramList[1]
    plotTypes = ["transmittance", "absorption coefficient", "cross section", "absorbance"]
    entryList = []
    for plotType in plotTypes:
        entryList.append(Entry(plotType, nextFunction=createPlot, functionParams=[obj, plotType, title]))
    choosePlotTypeMenu = Menu('Choose plot type', entryList)
    choosePlotTypeMenu.displayMenu()
    return


def menuChooseLayerToPlot(empty=None):
    entryList = []
    for layer in genericAtmosphere:
        nextEntry = Entry(layer.name, nextFunction=menuChoosePlotType, functionParams=([[layer], layer.title]))
        entryList.append(nextEntry)
        nextEntry = Entry('%s and components' % layer.name, nextFunction=createObjAndComponents, functionParams=layer)
        entryList.append(nextEntry)
    plotLayerMenu = Menu('Plot layer', entryList)
    plotLayerMenu.displayMenu()
    return


def createObjAndComponents(obj):
    plotList = [obj]
    for item in obj:
        plotList.append(item)
    menuChoosePlotType([plotList, obj.title])


def menuEditComposition(layer):
    moleculeList = layer.returnMoleculeObjects()
    entryList = []
    for molecule in moleculeList:
        newEntry = Entry('%s : %s' % (molecule.name, molecule.concText),
                         functionParams=molecule, nextFunction=inputMoleculeComposition)
        entryList.append(newEntry)
    entryList.append(Entry('Add a new molecule(s)', nextFunction=createMolecule, functionParams=layer))
    editCompMenu = Menu('Choose a molecule to edit', entryList, previousMenu=menuEditParamsOrComp)
    editCompMenu.displayMenu()
    return


def menuEditParamsOrComp(layer):
    entryList = []
    editLayerParamsEntry = Entry('Edit layer parameters', nextFunction=menuEditLayerParam, functionParams=layer)
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


def receiveInput(inputText, validInputFunction, default=None):
    validInput = False
    while not validInput:
        userInput = input(inputText)
        if userInput == '':
            return validInputFunction(str(default))
        validInput = validInputFunction(userInput)
    return validInputFunction(userInput)


def validMoleculeName(userInput):
    if not userInput:
        return False
    if userInput.strip().lower() == 'help':
        util.displayAllMolecules()
        return False
    elif userInput in pyrad.MOLECULE_ID:
        return userInput
    else:
        print('Invalid molecule name. %s' % (util.underlineMagenta('Please try again.')))
        return False


def validPressure(userInput):
    if not userInput:
        return False
    try:
        splitInput = validValueAndUnits.match(userInput)
    except AttributeError:
        print('Invalid input for pressure. Example: %s. %s'
              % (util.limeText('1.35atm'), util.underlineMagenta('Please try again.')))
        return False
    unit = splitInput.group(5)
    if not unit:
        unit = 'mbar'
    unit = unit.lower()
    if unit not in PRESSURE_UNITS:
        print('Invalid units. Accepted units are %s.' % ', '.join(PRESSURE_UNITS))
        return False
    textNumber = ''
    for i in range(1, 5):
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
        print('Invalid input for concentration. Example: %s. %s'
              % (util.limeText('15ppb'), util.underlineMagenta('Please try again.')))
        return False
    unit = splitInput.group(5)
    if not unit:
        unit = 'concentration'
    if unit not in COMPOSITION_UNITS:
        print('Invalid units. Accepted units are %s.'
              % (util.limeText(', '.join(COMPOSITION_UNITS))))
        return False
    textNumber = ''
    for i in range(1, 5):
        if splitInput.group(i):
            textNumber += splitInput.group(i)
    value = float(textNumber)
    if value <= 0:
        print('Concentration must be greater than 0')
        return False
    return value, unit


def validTemperature(userInput):
    if not userInput:
        return False
    try:
        splitInput = validValueAndUnits.match(userInput)
    except AttributeError:
        print('Invalid input for temperature. Example: %s. %s'
              % (util.limeText('20C'), util.underlineMagenta('Please try again.')))
        return False
    unit = splitInput.group(5)
    if not unit:
        unit = 'K'
    unit = unit.upper()[0]
    if unit not in TEMPERATURE_UNITS:
        print('Invalid units. Accepted units are %s.'
              % (util.limeText(', '.join(TEMPERATURE_UNITS))))
        return False
    textNumber = ''
    for i in range(1, 5):
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
        print('Invalid input for range. Example: %s. %s'
              % (util.limeText('150cm-1'), util.underlineMagenta('Please try again.')))
        return False
    unit = splitInput.group(5)
    if not unit:
        unit = 'cm-1'
    unit = unit.lower()
    if unit not in RANGE_UNITS:
        print('Invalid units. Accepted units are %s'
              % (util.limeText(', '.join(RANGE_UNITS))))
        return False
    textNumber = ''
    for i in range(1, 5):
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
        print('Invalid input for depth. Example: %s. %s'
              % (util.limeText('10cm'), util.underlineMagenta('Please try again.')))
        return False
    unit = splitInput.group(5)
    if not unit:
        unit = 'cm'
    unit = unit.lower()
    if unit not in DEPTH_UNITS:
        print('Invalid units. Accepted units are %s'
              % (util.limeText(', '.join(DEPTH_UNITS))))
        return False
    textNumber = ''
    for i in range(1, 5):
        if splitInput.group(i):
            textNumber += splitInput.group(i)
    value = float(textNumber)
    return value, unit


DEPTH_UNITS = ['cm', 'in', 'inches', 'ft', 'feet', 'meter', 'm']
PRESSURE_UNITS = ['atm', 'bar', 'mbar', 'pa']
TEMPERATURE_UNITS = ['K', 'C', 'F']
RANGE_UNITS = ['um', 'cm-1']
COMPOSITION_UNITS = ['ppm', 'ppb', '%', 'percentage', 'perc', 'concentration']

while True:
    menuMain()
