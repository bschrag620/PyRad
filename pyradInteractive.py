###code.interact(local=dict(globals(), **locals()))
import code
import pyrad
import pyradUtilities as util
import re
import pyradPlanck

existingAtmosphere = False
validValueAndUnits = re.compile('([-])?(\d+)?([.])?(\d+)?(\S+)?')
genericAtmosphere = pyrad.Atmosphere('holding atm for pyrad interactive')


class Menu:
    def __init__(self, title, entries, previousMenu=None, menuParams=None, hint=None):
        self.title = title
        self.entries = entries
        self.previousMenu = previousMenu
        self.menuParams = menuParams
        self.hint = hint

    def displayMenu(self):
        titleStr = '\t' + self.title
        while len(titleStr) < 60:
            titleStr += ' '
        print('\n%s' % util.underlineCyan(titleStr))
        
        i = 1
        validEntry = ['x']
        for entry in self.entries:
            if entry.selectionKey:
                validEntry.append(entry.selectionKey.lower())
                print(' %s)   %s' % (util.magentaText(entry.selectionKey.upper()), entry.name))
            else:
                validEntry.append(str(i))
                print(' %s)   %s' % (util.magentaText(i), entry.name))
            i += 1
        if 'Main' not in self.title:
            print(' %s   Previous menu' % util.magentaText('B)'))
            validEntry.append('b')
        print(' %s   Exit' % util.magentaText('X)'))
        validChoice = False
        if self.hint:
            print(util.limeText('**' + self.hint))

        while not validChoice:
            userInput = input('Choose an option: ')
            if userInput.lower() == 'x':
                print('Goodbye')
                exit(1)
            elif userInput.lower() == 'b' and 'Main' not in self.title:
                return
            elif userInput in validEntry:
                try:
                    userChoice = self.entries[int(userInput) - 1]
                except ValueError:
                    # this means the user entered a letter but it is in the valid choices
                    # roll through the entries to find the matching choice
                    userChoice = next(filter(lambda x: x.selectionKey and x.selectionKey.lower() == userInput.lower(), self.entries))
                if userChoice.nextFunction:
                    userChoice.nextFunction(userChoice.functionParams)
                    validChoice = True
                elif userChoice.nextMenu:
                    userChoice.nextMenu.displayMenu()
                    validChoice = True
            else:
                print('Invalid entry. Try again.')


class Entry:
    def __init__(self, text, nextMenu=None, nextFunction=None, functionParams=None, previousMenu=None, selectionKey=False):
        self.name = text
        self.nextMenu = nextMenu
        self.nextFunction = nextFunction
        self.functionParams = functionParams
        self.previousMenu = previousMenu
        self.selectionKey = selectionKey


def createPlot(params):
    plotList = params['plots']
    plotType = params['plotType']
    plotTitle = params['title']
    pyrad.plot(plotType, plotTitle, plotList)
    return


def createLayer(atmosphere):
    defaultLayerName = atmosphere.nextLayerName()
    getLayerName = '\n%s\n' \
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
        if moleculeName in XSC_LIST:
            params = {'layer': layer,
                        'xsc': moleculeName}
            selectXscFile(params)
        else:
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


def addXscToLayer(params):
    layer = params['layer']
    xsc = params['xsc']
    file = params['file']
    concentration, units = inputMoleculeComposition()
    tempdict = {units: concentration}
    layer.addMolecule({xsc: file}, **tempdict)
    return

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
    depth = receiveInput(getDepth, validDepth, default=default)
    return depth


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
    temperature = receiveInput(getTemperature, validTemperature, default=default)
    return temperature


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
    pressure = receiveInput(getPressure, validPressure, default=default)
    return pressure


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
        rangeMin = receiveInput(getRangeMin, validRange, default=defaultMin)
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
        rangeMax = receiveInput(getRangeMax, validRange, default=defaultMax)
        if rangeMax <= rangeMin:
            print('Range min must be %s than range min of %s' % (util.magentaText('greater'), util.cyanText(rangeMin)))
    return rangeMin, rangeMax


def validNumber(userInput):
    try:
        return float(userInput)
    except ValueError:
        return False


def inputPlanckRange(units):
    rangeMin = -1
    rangeMax = -1
    text = '%s\nUnits are %s. Scientific notation is accepted (1e14):' \
           % (util.underlineCyan('Enter the minimum range of the planck spectrum.'),util.limeText(units))
    while rangeMin < 0:
        rangeMin = receiveInput(text, validNumber)
        if rangeMin < 0:
            print('Range min must be %s than zero' % util.magentaText('greater'))
    text = '%s\nUnits are %s. Scientific notation is accepted (1e14):' \
           % (util.underlineCyan('Enter the maximum range of the planck spectrum.'), util.limeText(units))
    while rangeMax <= rangeMin:
        rangeMax = receiveInput(text, validNumber)
        if rangeMax <= rangeMin:
            print('Range min must be %s than range min of %s' % (util.magentaText('greater'), util.cyanText(rangeMin)))
    return rangeMin, rangeMax


def inputMoleculeName(default=None):
    if not default:
        default = 'co2'
    text = 'Enter the short molecule name.\t\t\t'
    moleculeName = receiveInput('%s\n'
                                'For a full list of options, type %s. For a list of cross-section only molecules, type %s. If no value given, %s will be used: '
                                % (util.underlineCyan(text),
                                   util.magentaText('help'),
                                   util.magentaText('xsc'),
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


def inputPlanckTemps():
    text = 'Enter the temperature of the planck curves.\t\t\t'
    tempList = receiveMultiInput('%s\n'
                            'If no units entered, temperature will be assumed %s.\n'
                            'Other valid units are %s .Multiple temperatures can be separated with a comma: '
                            % (util.underlineCyan(text),
                               util.limeText('K'),
                               util.limeText('C and F')), validTemperature)
    return tempList


def menuChooseLayerToEdit(empty=None):
    entryList = []
    for layer in genericAtmosphere:
        nextEntry = Entry(layer.name, nextFunction=menuEditParamsOrComp, functionParams=layer)
        entryList.append(nextEntry)
    editLayerMenu = Menu('Edit layer', entryList)
    editLayerMenu.displayMenu()
    return


def menuChooseTransmission(plotType):
    entryList = []
    for layer in genericAtmosphere:
        params = {'plots': [layer], 'plotType': plotType, 'title': layer.title}
        nextEntry = Entry(layer.name, nextFunction=createTransmission, functionParams=params)
        entryList.append(nextEntry)
        plotList, title = createObjAndComponents(layer)
        params = {'plots': plotList, 'plotType': plotType, 'title': title}
        nextEntry = Entry('%s and components' % layer.name, nextFunction=createTransmission, functionParams=params)
        entryList.append(nextEntry)
    transmissionMenu = Menu('Choose which layers to plot transmission', entryList)
    transmissionMenu.displayMenu()
    return


def createPlanckCurves(plotType):
    plotList = inputPlanckTemps()
    rangeMin, rangeMax = inputPlanckRange(plotType)
    pyrad.plotSpectrum(title='Planck spectrums', rangeMin=rangeMin, rangeMax=rangeMax, planckTemperatureList=plotList, planckType=plotType)
    return


def menuPlanckType(empty=None):
    entryList = []
    wavenumber = 'wavenumber'
    hertz = 'Hz'
    wavelength = 'wavelength'
    entryList.append(Entry('By %s (cm-1)' % wavenumber, nextFunction=createPlanckCurves, functionParams=wavenumber))
    entryList.append(Entry('By %s (um)' % wavelength, nextFunction=createPlanckCurves, functionParams=wavelength))
    entryList.append(Entry('By %s (s-1)' % hertz, nextFunction=createPlanckCurves, functionParams=hertz))
    planckTypeMenu = Menu('Choose planck type', entryList)
    planckTypeMenu.displayMenu()
    return


def createTransmission(params):
    layer = params['plots'][0]
    text = 'A plot for transmission requires an initial surface temperature.\n'\
           'Please choose a temperature different from the layer temperature of %sK:' % util.limeText(layer.T)
    temperature = receiveMultiInput(text, validTemperature)
    objList = []
    for item in params['plots']:
        objList.append(item)
    temperature.append(layer.T)
    pyrad.plotSpectrum(layer, objList=objList,
                       surfaceSpectrum=pyradPlanck.planckWavenumber(layer.xAxis, temperature[0]),
                       planckTemperatureList=temperature)
    return


def menuChoosePlotType(empty=None):
    entryList = []
    entryList.append(Entry('transmittance', nextFunction=menuChooseLayerToPlot, functionParams='transmittance'))
    entryList.append(Entry('absorption coefficient', nextFunction=menuChooseLayerToPlot, functionParams='absorption coefficient'))
    entryList.append(Entry('cross section', nextFunction=menuChooseLayerToPlot, functionParams='cross section'))
    entryList.append(Entry('absorbance', nextFunction=menuChooseLayerToPlot, functionParams='absorbance'))
    entryList.append(Entry('optical depth', nextFunction=menuChooseLayerToPlot, functionParams='optical depth'))
    entryList.append(Entry('line survey', nextFunction=menuChooseLayerToPlot, functionParams='line survey'))
    entryList.append(Entry('transmission', nextFunction=menuChooseTransmission, functionParams='transmission'))
    choosePlotTypeMenu = Menu('Choose plot type', entryList)
    choosePlotTypeMenu.displayMenu()
    return


def menuChooseLayerToPlot(plotType):
    entryList = []
    for layer in genericAtmosphere:
        params = {'plots': [layer], 'plotType': plotType, 'title': layer.title}
        nextEntry = Entry(layer.name, nextFunction=createPlot, functionParams=params)
        entryList.append(nextEntry)
        plotList, title = createObjAndComponents(layer)
        params = {'plots': plotList, 'plotType': plotType, 'title': title}
        nextEntry = Entry('%s and components' % layer.name, nextFunction=createPlot, functionParams=params)
        entryList.append(nextEntry)
    plotLayerMenu = Menu('Plot layer', entryList)
    plotLayerMenu.displayMenu()
    return


def createObjAndComponents(obj):
    plotList = [obj]
    for item in obj:
        plotList.append(item)
    return plotList, obj.title


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
    entries = []
    entries.append(Entry("Create new gas cell", nextFunction=createLayer, functionParams=genericAtmosphere))
    entries.append(Entry("Edit/duplicate gas cell", nextFunction=menuChooseLayerToEdit))
    entries.append(Entry("Plot gas cell", nextFunction=menuChoosePlotType))
    entries.append(Entry('Plot planck curves', nextFunction=menuPlanckType))
    mainMenu = Menu('Main menu', entries)
    mainMenu.displayMenu()
    return


def menuAddXscToLayer(atm):
    entries = []
    for layer in atm:
        entries.append(Entry(layer.name, nextFunction=menuSelectXscMolecule, functionParams=layer))
    menu = Menu('Choose layer to add cross-section to', entries)
    menu.displayMenu()
    return


def menuSelectXscMolecule(layer):
    entries = []
    for xsc in XSC_LIST:
        entries.append(Entry(xsc, nextFunction=selectXscFile, functionParams={'layer': layer, 'xsc': xsc}))
    entries.append(Entry('Manually enter molecule name based on HITRAN list', nextFunction=inputXscName, functionParams=layer))
    menu = Menu('Choose xsc from list', entries, hint='Full list is available at https://hitran.org/data/suppl/xsec/cross_section_data/')
    menu.displayMenu()
    return


def selectXscFile(params):
    layer = params['layer']
    xsc = params['xsc']

    if 'sort' not in params:
        params['sort'] = 'TEMP'
    sort = params['sort']
    entries = []
    unsortedFiles = util.returnXscFilesInDirectory(xsc)
    if unsortedFiles is False or unsortedFiles == []:
        question = "Either that directory doesn't exist or it was empty. Would you like to try downloading the data from HITRAN?"
        if receiveInput(question, validYorN) == 'y':
            filepath = util.downloadXscZipFile(xsc)
            util.unzipFile(filepath)
            util.mergeXsc(xsc)
            selectXscFile(params)
        else:
            return
    else:
        unsortedValues = list(map(lambda file: util.parseXscFileName(file), unsortedFiles))
        if sort == 'TEMP':
            sortedValues = sorted(unsortedValues, key = lambda i: (float(i['TEMP']), float(i['PRESSURE'])))
        elif sort == 'PRESSURE':
            sortedValues = sorted(unsortedValues, key = lambda i: (float(i['PRESSURE']), float(i['TEMP'])))
        for v in sortedValues:
            if sort == 'TEMP':
                displayName = 'Temp: %s  -- Pressure: %s  --  Range: %s' % (util.limeText(v['TEMP'] + 'K'), util.cyanText(v['PRESSURE'] + 'Torr'), util.magentaText(v['RANGE'] + 'cm-1'))
            elif sort == 'PRESSURE':
                displayName = 'Pressure: %s  -- Temp: %s  --  Range: %s' % (util.cyanText(v['PRESSURE'] + 'Torr'), util.limeText(v['TEMP'] + 'K'), util.magentaText(v['RANGE'] + 'cm-1'))
            entries.append(Entry(displayName, nextFunction=addXscToLayer, functionParams={'layer': layer, 'file': v['LONG_FILENAME'], 'xsc': xsc}))
        pressureParams = params.copy()
        pressureParams.update({'sort': 'PRESSURE'})
        entries.append(Entry('Sort by pressure', nextFunction=selectXscFile, functionParams=pressureParams, selectionKey='P'))
        tempParams = params.copy()
        tempParams.update({'sort': 'TEMP'})
        entries.append(Entry('Sort by temperature', nextFunction=selectXscFile, functionParams=tempParams, selectionKey='T'))
        menu = Menu('Choose file to use', entries, hint='Layer P and T will be adjusted according to the xsc file')
        menu.displayMenu()
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
    while validInput is False:
        userInput = input('\n%s' % inputText)
        if userInput == '':
            return validInputFunction(str(default))
        validInput = validInputFunction(userInput)
    return validInput


def receiveMultiInput(inputText, validInputFunction, default=None):
    validInput = False
    while not validInput:
        userInput = input(inputText)
        inputList = userInput.replace(' ', '').split(',')
        testInput = True
        for item in inputList:
            if not validInputFunction(item):
                print('%s not recognized as valid. Please try again.')
                testInput = False
        if testInput:
            return inputList


def validYorN(userInput):
    if not userInput:
        return False
    if userInput.strip().lower()[0] == 'y' or userInput.strip().lower()[0] == 'n':
        return userInput
    else:
        print('Invalid option. Please type "%s" or "%s"' % (util.underlineMagenta('y'), util.underlineMagenta('n')))
        return False


def validMoleculeName(userInput):
    if not userInput:
        return False
    if userInput.strip().lower() == 'help':
        util.displayAllMolecules()
        return False
    elif userInput.strip().lower() == 'xsc':
        print(util.underlineMagenta('Punctuation matters...'))
        print(', '.join(XSC_LIST))
        return False
    elif userInput in pyrad.MOLECULE_ID or userInput in XSC_LIST:
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
    return pyrad.convertPressure(value, unit)


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
    return pyrad.convertTemperature(value, unit)


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
    return pyrad.convertRange(value, unit)


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
    return pyrad.convertLength(value, unit)


DEPTH_UNITS = ['cm', 'in', 'inches', 'ft', 'feet', 'meter', 'm']
PRESSURE_UNITS = ['atm', 'bar', 'mbar', 'pa']
TEMPERATURE_UNITS = ['K', 'C', 'F']
RANGE_UNITS = ['um', 'cm-1']
COMPOSITION_UNITS = ['ppm', 'ppb', '%', 'percentage', 'perc', 'concentration']
XSC_LIST = ['CFC-11', 'CFC-12', 'CFC-13', 'CFC-113', 'CFC-113a', 'CFC-114', 'CFC-114a', 'CFC-115',
            'HCFC-21', 'HCFC-22', 'HCFC-123', 'HCFC-123a', 'HCFC-124', 'HCFC-141b', 'HCFC-142b', 'HCFC-225ca', 'HCFC-225cb',
            'HFC-32', 'HFC-125', 'HFC-134', 'HFC-134a', 'HFC-143a', 'HFC-152a', 'HFE-356mff2']

while True:
    menuMain()