import pyradClasses
import pyradUtilities as util
import re
import pyradPlanck

existingAtmosphere = False
validValueAndUnits = re.compile('([-])?(\d+)?([.])?(\d+)?(\S+)?')
genericAtmosphere = pyradClasses.Atmosphere('holding atm for pyrad interactive')
pyradClasses.settings.changeSetting('hi')


class Menu:
    def __init__(self, title, entries, previousMenu=None, menuParams=None):
        self.title = title
        self.entries = entries
        self.previousMenu = previousMenu
        self.menuParams = menuParams

    def displayMenu(self):
        titleStr = '\t%s\tdetail: %s' % (self.title, pyradClasses.settings.setting)
        while len(titleStr) < 60:
            titleStr += ' '
        print('\n%s' % util.underlineCyan(titleStr))
        i = 1
        validEntry = ['x']
        for entry in self.entries:
            validEntry.append(str(i))
            print(' %s)   %s' % (util.magentaText(i), entry.name))
            i += 1
        if 'Main' not in self.title:
            print(' %s   Previous menu' % util.magentaText('B)'))
            validEntry.append('b')
        if 'settings' not in self.title:
            print(' %s   Settings' % util.magentaText('S)'))
            validEntry.append('s')
        print(' %s   Exit' % util.magentaText('X)'))
        validChoice = False
        while not validChoice:
            userInput = input('Choose an option: ')
            if userInput.lower() == 'x':
                print('Goodbye')
                exit(1)
            elif userInput.lower() == 'b' and 'main' not in self.title.lower():
                self.previousMenu()
            elif userInput.lower() == 's' and 'settings' not in self.title.lower():
                settingsMenu(self)
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

    def displayMultiChoiceMenu(self):
        titleStr = '\t%s\tdetail: %s' % (self.title, pyradClasses.settings.setting)
        while len(titleStr) < 60:
            titleStr += ' '
        print('\n%s' % util.underlineCyan(titleStr))
        i = 1
        validEntry = ['x']
        for entry in self.entries:
            validEntry.append(str(i))
            print(' %s)   %s' % (util.magentaText(i), entry.name))
            i += 1
        if 'Main' not in self.title:
            print(' %s   Previous menu' % util.magentaText('B)'))
            validEntry.append('b')
        if 'settings' not in self.title:
            print(' %s   Settings' % util.magentaText('S)'))
            validEntry.append('s')
        print(' %s   Exit' % util.magentaText('X)'))
        validChoice = False
        while not validChoice:
            userInput = input('Choose an option: ')
            if userInput.lower() == 'x':
                print('Goodbye')
                exit(1)
            elif userInput.lower() == 'b' and 'Main' not in self.title:
                return
            elif userInput.lower() == 's' and 'settings' not in self.title.lower():
                settingsMenu(self)
            else:
                inputs = userInput.split(',')
                allValid = True
                userChoices = []
                for i in inputs:
                    print('i: %s' % i)
                    if i.strip() not in validEntry:
                        allValid = False
                    else:
                        userChoices.append(self.entries[int(i) - 1])
                if allValid:
                    nextFunction = userChoices[0].nextFunction
                    nextFunction(userChoices)
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
    plotList = params['plots']
    plotType = params['plotType']
    plotTitle = params['title']
    pyradClasses.plot(plotType, plotTitle, plotList)
    menuChoosePlotType()


def plotPlanetSpectrum(values):
    pList = []
    for p in values['profiles']:
        planet = pyradClasses.loadEmptyPlanet(p.name)
        pList.append(planet)
    if values['height'] == -2.71828:
        pyradClasses.plotPlanetSpectrum(pList, direction=values['direction'], verify=False)
    else:
        pyradClasses.plotPlanetSpectrum(pList, direction=values['direction'], height=values['height'], verify=False)
    chooseAtmTransferBuildProfile()


def plotPlanetSpectrumComponents(values):
    planet = pyradClasses.createCustomPlanet(values['profiles'].name)
    if values['height'] == -2.71828:
        pyradClasses.plotPlanetAndComponents(planet, direction=values['direction'], verify=False)
    else:
        pyradClasses.plotPlanetAndComponents(planet, direction=values['direction'], height=values['height'], verify=False)
    chooseAtmTransferBuildProfile()


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
        concentration, units = inputMoleculeComposition()
        tempdict = {units: concentration}
        molecule = layer.addMolecule(moleculeName, **tempdict)
        while pyradClasses.totalConcentration(layer) > 1:
            print("%s total concentration exceeds 100%%" % util.magentaText('***\tWARNING\t***'))
            menuEditComposition(layer)
        validInput = False
        while not validInput:
            ask = input("Add another molecule to the layer %s " % util.magentaText('(y/n) :'))
            if ask.strip().lower() == 'y':
                validInput = True
            elif ask.strip().lower() == 'n':
                gasCellMenu()


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


def inputPlanckTemps():
    text = 'Enter the temperature of the planck curves.\t\t\t'
    tempList = receiveMultiInput('%s\n'
                            'If no units entered, temperature will be assumed %s.\n'
                            'Other valid units are %s .Multiple temperatures can be separated with a comma: '
                            % (util.underlineCyan(text),
                               util.limeText('K'),
                               util.limeText('C and F')), validTemperature)
    return tempList


def inputHeight(values):
    text = 'Enter the height to view transmission from.\t\t\t'
    height = receiveInput('%s\n'
                          'Units should be in %s. If no value entered, maximum atm height will be used: ' % (util.underlineCyan(text), util.limeText('km')), validNumber, default=-2.71828)
    values['height'] = height
    plotPlanetSpectrum(values)
    return


def inputHeightComponents(values):
    text = 'Enter the height to view transmission from.\t\t\t'
    height = receiveInput('%s\n'
                          'Units should be in %s. If no value entered, maximum atm height will be used: ' % (util.underlineCyan(text), util.limeText('km')), validNumber, default=-2.71828)
    values['height'] = height
    plotPlanetSpectrumComponents(values)
    return


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
    pyradClasses.plotSpectrum(title='Planck spectrums', rangeMin=rangeMin, rangeMax=rangeMax, planckTemperatureList=plotList, planckType=plotType)
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
    pyradClasses.plotSpectrum(layer, objList=objList,
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
    choosePlotTypeMenu = Menu('Choose plot type', entryList, previousMenu=gasCellMenu)
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
    plotLayerMenu = Menu('Plot layer', entryList, previousMenu=menuChoosePlotType)
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


def gasCellMenu(param=None):
    createLayerEntry = Entry("Create new gas cell", nextFunction=createLayer, functionParams=genericAtmosphere)
    editLayerEntry = Entry("Edit/duplicate gas cell", nextFunction=menuChooseLayerToEdit)
    plotLayerEntry = Entry("Plot gas cell", nextFunction=menuChoosePlotType)
    planckPlotEntry = (Entry('Plot planck curves', nextFunction=menuPlanckType))
    menuGasCell = Menu('Gas cell simulator', [createLayerEntry, editLayerEntry, plotLayerEntry, planckPlotEntry], previousMenu=menuMain)
    menuGasCell.displayMenu()
    return


def chooseDirection(profileList):
    lookUpEntry = Entry('Looking up', nextFunction=inputHeight, functionParams={'profiles': profileList,
                                                                                'direction': 'up'})
    lookDownEntry = Entry('Looking down', nextFunction=inputHeight, functionParams={'profiles': profileList,
                                                                                'direction': 'down'})
    menuChooseDirection = Menu('Choose direction to look', [lookUpEntry, lookDownEntry], previousMenu=plotAtmTransferMenu)
    menuChooseDirection.displayMenu()
    return


def chooseDirectionComponents(profile):
    lookUpEntry = Entry('Looking up', nextFunction=inputHeightComponents, functionParams={'profiles': profile,
                                                                                'direction': 'up'})
    lookDownEntry = Entry('Looking down', nextFunction=inputHeightComponents, functionParams={'profiles': profile,
                                                                                    'direction': 'down'})
    menuChooseDirection = Menu('Choose direction to look', [lookUpEntry, lookDownEntry],
                               previousMenu=plotProfileComponentsMenu)
    menuChooseDirection.displayMenu()
    return


def plotProfileComponentsMenu(param=None):
    entryList = []
    profileList = util.getProfileList()
    for profile in profileList:
        if util.profileComplete('%s %s' % (profile, pyradClasses.settings.setting)) and \
                util.molSpecProfile(profile, pyradClasses.settings.setting):
            entryList.append(Entry('%s' % profile[:-4], nextFunction=chooseDirectionComponents, functionParams=profile))
    menuAtmTransfer = Menu('Choose atmosphere', entryList, previousMenu=chooseAtmTransferBuildProfile)
    menuAtmTransfer.displayMenu()
    return


def buildProfile(profileList):
    for profile in profileList:
        print('Building %s on setting %s' % (profile.name, pyradClasses.settings.setting))
        planet = pyradClasses.createCustomPlanet(profile.name)
        moleculeSpecific = pyradClasses.yesOrNo("Store data by individual molecule? Doesn't require more time, "
                                                "just hard drive space %s:" % util.limeText('(y/n)'))
        overwrite = True
        if util.profileComplete(planet.folderPath):
            overwrite = pyradClasses.yesOrNo("Data for this profile and setting seems to exist.\n"
                                             "Do you wish to overwrite it? %s" % util.limeText('(y/n)'))
        if overwrite:
            planet.processLayers(verify=False, moleculeSpecific=moleculeSpecific)
    chooseAtmTransferBuildProfile()


def buildProfilesMenu(param=None):
    entryList = []
    profileList = util.getProfileList()
    for profile in profileList:
        entryList.append(Entry('%s' % profile[:-4], nextFunction=buildProfile, functionParams=profile))
    menuBuildProfile = Menu('Which profile(s) do you want to build', entryList, previousMenu=chooseAtmTransferBuildProfile)
    menuBuildProfile.displayMultiChoiceMenu()
    return


def chooseAtmTransferBuildProfile(param=None):
    plotProfilesEntry = Entry('Plot profile(s)', nextFunction=plotAtmTransferMenu)
    plotProfileAndComponents = Entry('Plot single atm and components', nextFunction=plotProfileComponentsMenu)
    buildProfilesEntry = Entry('Build profile(s)', nextFunction=buildProfilesMenu)
    menuChooseProfileBuild = Menu('Plot or build profiles', [plotProfilesEntry, plotProfileAndComponents, buildProfilesEntry], previousMenu=menuMain)
    menuChooseProfileBuild.displayMenu()
    return


def plotAtmTransferMenu(param=None):
    profileList = util.getCompletedProfileList(pyradClasses.settings.setting)
    entryList = []
    for profile in profileList:
        entryList.append(Entry('%s' % profile, nextFunction=chooseDirection, functionParams=profile))
    menuAtmTransfer = Menu('Choose atmosphere', entryList, previousMenu=menuMain)
    menuAtmTransfer.displayMultiChoiceMenu()
    return


def settingsMenu(previousMenu):
    lowSetting = Entry('low (intensity > 1E-21)', nextFunction=pyradClasses.settings.changeSetting, functionParams='low')
    midSetting = Entry('mid (intensity > 1E-28)', nextFunction=pyradClasses.settings.changeSetting, functionParams='mid')
    hiSetting = Entry('hi (all absorption lines)', nextFunction=pyradClasses.settings.changeSetting, functionParams='hi')
    menuSettings = Menu('Choose level of detail', [lowSetting, midSetting, hiSetting], previousMenu=menuMain)
    return menuSettings


def menuMain():
    gasCellEntry = Entry('Gas cell simulator', nextFunction=gasCellMenu)
    atmosphereTransferEntry = Entry('Atmosphere transmission', nextFunction=chooseAtmTransferBuildProfile)
    mainMenu = Menu('Main menu', [gasCellEntry, atmosphereTransferEntry])
    return mainMenu



def duplicateObj(obj):
    newObj = obj.returnCopy()
    if isinstance(newObj, pyradClasses.Layer):
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


def validMoleculeName(userInput):
    if not userInput:
        return False
    if userInput.strip().lower() == 'help':
        util.displayAllMolecules()
        return False
    elif userInput in pyradClasses.MOLECULE_ID:
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
    return pyradClasses.convertPressure(value, unit)


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
    return pyradClasses.convertTemperature(value, unit)


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
    return pyradClasses.convertRange(value, unit)


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
    return pyradClasses.convertLength(value, unit)


DEPTH_UNITS = ['cm', 'in', 'inches', 'ft', 'feet', 'meter', 'm']
PRESSURE_UNITS = ['atm', 'bar', 'mbar', 'pa']
TEMPERATURE_UNITS = ['K', 'C', 'F']
RANGE_UNITS = ['um', 'cm-1']
COMPOSITION_UNITS = ['ppm', 'ppb', '%', 'percentage', 'perc', 'concentration']

menu = menuMain()
while True:
    menu = menu.displayMenu()