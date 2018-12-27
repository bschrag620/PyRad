import pyradClasses
import pyradUtilities as util
import re
import pyradPlanck

existingAtmosphere = False
validValueAndUnits = re.compile('([-])?(\d+)?([.])?(\d+)?(\S+)?')
genericAtmosphere = pyradClasses.Atmosphere('holding atm for pyrad interactive')
pyradClasses.settings.changeSetting('hi')


class Menu:
    def __init__(self, title, entries, previousMenu=None, menuParams=None, multiChoice=False):
        self.title = title
        self.entries = entries
        self.previousMenu = previousMenu
        self.menuParams = menuParams
        self.multiChoice = multiChoice

    def displayMenu(self):
        titleStr = '   %s   detail: %s' % (self.title, pyradClasses.settings.setting)
        while len(titleStr) < 60:
            titleStr += ' '
        print('\n%s' % util.underlineCyan(titleStr))
        i = 1
        validEntry = ['x']
        for entry in self.entries:
            validEntry.append(str(i))
            print(' %s)   %s' % (util.magentaText(i), entry.name))
            i += 1
        text = ' %s)   Exit' % util.magentaText('X')
        if 'Main' not in self.title:
            text += '\t%s)  Back' % util.magentaText('B')
            validEntry.append('b')
        if self.title != 'Choose level of detail':
            text += '\t%s)  Settings' % util.magentaText('S')
            validEntry.append('s')
        print(text)
        validChoice = False
        while not validChoice:
            userInput = input('Choose an option: ')
            if userInput.lower() == 'x':
                print('Goodbye')
                exit(1)
            elif userInput.lower() == 'b' and 'main' not in self.title.lower():
                if type(self.previousMenu) == Menu:
                    return self.previousMenu
                return self.previousMenu()
            elif userInput.lower() == 's' and 'settings' not in self.title.lower():
                return settingsMenu(previousMenu=self)
            elif not self.multiChoice and userInput in validEntry:
                userChoice = self.entries[int(userInput) - 1]
                if userChoice.nextFunction:
                    return userChoice.nextFunction(userChoice.functionParams)
                elif userChoice.nextMenu:
                    return userChoice.nextMenu(userChoice.functionParams)
            elif self.multiChoice:
                inputs = userInput.split(',')
                allValid = True
                userChoices = []
                for i in inputs:
                    if i.strip() not in validEntry:
                        allValid = False
                    else:
                        userChoices.append(self.entries[int(i) - 1])
                if allValid:
                    if userChoices[0].nextFunction:
                        nextFunction = userChoices[0].nextFunction
                        return nextFunction(userChoices)
                    else:
                        nextMenu = userChoices[0].nextMenu
                        tempMenu = nextMenu(userChoices)
                        return tempMenu
                else:
                    print('Invalid entry. Try again.')
            else:
                print('Invalid entry. Try again.')

    def displayMultiChoiceMenu(self):
        titleStr = '%sdetail: %s' % (self.title, pyradClasses.settings.setting)
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
            userInput = input('Choose an option(s). Separate multiple choices with a comma: ')
            if userInput.lower() == 'x':
                print('Goodbye')
                exit(1)
            elif userInput.lower() == 'b' and 'main' not in self.title.lower():
                return self.previousMenu()
            elif userInput.lower() == 's' and 'settings' not in self.title.lower():
                return settingsMenu(previousMenu=self)
            else:
                inputs = userInput.split(',')
                allValid = True
                userChoices = []
                for i in inputs:

                    if i.strip() not in validEntry:
                        allValid = False
                    else:
                        userChoices.append(self.entries[int(i) - 1])
                if allValid:
                    if userChoices[0].nextFunction:
                        nextFunction = userChoices[0].nextFunction
                        return nextFunction(userChoices)
                    else:
                        nextMenu = userChoices[0].nextMenu
                        tempMenu = nextMenu(userChoices)
                        return tempMenu
                else:
                    print('Invalid entry. Try again.')


class Entry:
    def __init__(self, text, nextMenu=None, nextFunction=None, functionParams=None, previousMenu=None):
        self.name = text
        self.nextMenu = nextMenu
        self.nextFunction = nextFunction
        self.functionParams = functionParams
        if previousMenu:
            self.previousMenu = previousMenu()


def changeSettings(params):
    value = params[0]
    previousMenu = params[1]
    pyradClasses.settings.changeSetting(value)
    return previousMenu


def loadTheme(params):
    value = params[0]
    previousMenu = params[1]
    pyradClasses.theme.loadTheme(value)
    return settingsMenu(previousMenu)


def createPlot(params):
    plotList = params['plots']
    plotType = params['plotType']
    plotTitle = params['title']
    pyradClasses.plot(plotType, plotTitle, plotList)
    return menuChoosePlotType()


def plotPlanetSpectrum(values):
    pList = []
    for p in values['profiles']:
        pList.append(p.name)
    if values['height'] == -2.71828:
        pyradClasses.plotPlanetSpectrum(pList, direction=values['direction'], verify=False)
    else:
        pyradClasses.plotPlanetSpectrum(pList, direction=values['direction'], height=values['height'], verify=False)
    return chooseAtmTransferBuildProfile()


def plotPlanetSpectrumComponents(values):
    #planet = pyradClasses.loadEmptyPlanet(values['profiles'][0].name)
    if values['height'] == -2.71828:
        pyradClasses.plotPlanetAndComponents(values['profiles'][0].name, direction=values['direction'], verify=False)
#        pyradClasses.plotPlanetAndComponents(values['profiles'], direction=values['direction'], verify=False)
    else:
        pyradClasses.plotPlanetAndComponents(values['profiles'][0].name, direction=values['direction'], height=values['height'], verify=False)
    return chooseAtmTransferBuildProfile()


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
    return createMolecule(layer)


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
                return gasCellMenu()


def menuEditLayerParam(layer):
    editDepth = Entry('Depth', nextFunction=editLayerDepth, functionParams=layer)
    editRange = Entry('Min or max range', nextFunction=editLayerRange, functionParams=layer)
    editTemperature = Entry('Temp   erature',nextFunction=editLayerTemperature, functionParams=layer)
    editPressure = Entry('Pressure', nextFunction=editLayerPressure, functionParams=layer)
    entryList = [editDepth, editRange, editTemperature, editPressure]
    menu = Menu('Choose the parameter to edit for %s' % layer.name, entryList,
                previousMenu=menuEditParamsOrComp, menuParams=layer)
    return menu


def editLayerDepth(layer):
    print('Current %s for %s is : %s\n'
          % (util.limeText('depth'), util.limeText(layer.name), util.cyanText('%scm' % layer.depth)))
    depth = inputLayerDepth(default=layer.depth)
    layer.changeDepth(depth)
    return menuEditLayerParam(layer)


def editLayerTemperature(layer):
    print('Current %s for %s is : %s\n'
          % (util.limeText('temperature'), util.limeText(layer.name), util.cyanText('%sK' % layer.T)))
    temperature = inputLayerTemperature(default=layer.T)
    layer.changeTemperature(temperature)
    return menuEditLayerParam(layer)


def editLayerPressure(layer):
    print('Current %s for %s is : %s\n'
          % (util.limeText('pressure'), util.limeText(layer.name), util.cyanText('%smbar' % layer.P)))
    pressure = inputLayerPressure(default=layer.P)
    layer.changePressure(pressure)
    return menuEditLayerParam(layer)


def editLayerRange(layer):
    print('Current %s for %s is %s\n'
          % (util.limeText('range'), util.limeText(layer.name),
             util.cyanText('%s-%scm-1' % (layer.rangeMin, layer.rangeMax))))
    rangeMin, rangeMax = inputLayerRange(defaultMin=layer.rangeMin, defaultMax=layer.rangeMax)
    layer.changeRange(rangeMin, rangeMax)
    return menuEditLayerParam(layer)


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
    values['height'] = height * 10**5
    return plotPlanetSpectrum(values)


def inputHeightComponents(values):
    text = 'Enter the height to view transmission from.\t\t\t'
    height = receiveInput('%s\n'
                          'Units should be in %s. If no value entered, maximum atm height will be used: ' % (util.underlineCyan(text), util.limeText('km')), validNumber, default=-2.71828)
    values['height'] = height * 10**5
    return plotPlanetSpectrumComponents(values)


def menuChooseLayerToEdit(empty=None):
    entryList = []
    for layer in genericAtmosphere:
        nextEntry = Entry(layer.name, nextFunction=menuEditParamsOrComp, functionParams=layer)
        entryList.append(nextEntry)
    editLayerMenu = Menu('Edit layer', entryList)
    return editLayerMenu


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
    return transmissionMenu


def createPlanckCurves(plotType):
    plotList = inputPlanckTemps()
    rangeMin, rangeMax = inputPlanckRange(plotType)
    pyradClasses.plotSpectrum(title='Planck spectrums', rangeMin=rangeMin, rangeMax=rangeMax, planckTemperatureList=plotList, planckType=plotType)
    return menuMain()


def menuPlanckType(empty=None):
    entryList = []
    wavenumber = 'wavenumber'
    hertz = 'Hz'
    wavelength = 'wavelength'
    entryList.append(Entry('By %s (cm-1)' % wavenumber, nextFunction=createPlanckCurves, functionParams=wavenumber))
    entryList.append(Entry('By %s (um)' % wavelength, nextFunction=createPlanckCurves, functionParams=wavelength))
    entryList.append(Entry('By %s (s-1)' % hertz, nextFunction=createPlanckCurves, functionParams=hertz))
    planckTypeMenu = Menu('Choose planck type', entryList, previousMenu=menuMain)
    return planckTypeMenu


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
    return menuChoosePlotType()


def menuChoosePlotType(empty=None):
    entryList = []
    entryList.append(Entry('transmittance', nextMenu=menuChooseLayerToPlot, functionParams='transmittance'))
    entryList.append(Entry('absorption coefficient', nextMenu=menuChooseLayerToPlot, functionParams='absorption coefficient'))
    entryList.append(Entry('cross section', nextMenu=menuChooseLayerToPlot, functionParams='cross section'))
    entryList.append(Entry('absorbance', nextMenu=menuChooseLayerToPlot, functionParams='absorbance'))
    entryList.append(Entry('optical depth', nextMenu=menuChooseLayerToPlot, functionParams='optical depth'))
    entryList.append(Entry('line survey', nextMenu=menuChooseLayerToPlot, functionParams='line survey'))
    entryList.append(Entry('transmission', nextMenu=menuChooseTransmission, functionParams='transmission'))
    choosePlotTypeMenu = Menu('Choose plot type', entryList, previousMenu=gasCellMenu)
    return choosePlotTypeMenu


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
    return plotLayerMenu


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
    return editCompMenu


def menuEditParamsOrComp(layer):
    entryList = []
    editLayerParamsEntry = Entry('Edit layer parameters', nextFunction=menuEditLayerParam, functionParams=layer)
    entryList.append(editLayerParamsEntry)
    duplicateLayerEntry = Entry('Duplicate layer', nextFunction=duplicateObj, functionParams=layer)
    entryList.append(duplicateLayerEntry)
    editCompositionEntry = Entry('Edit composition', nextFunction=menuEditComposition, functionParams=layer)
    entryList.append(editCompositionEntry)
    chooseParamsMenu = Menu('Edit or duplicate', entryList, previousMenu=menuChooseLayerToEdit)
    return chooseParamsMenu


def gasCellMenu(param=None):
    createLayerEntry = Entry("Create new gas cell", nextFunction=createLayer, functionParams=genericAtmosphere)
    editLayerEntry = Entry("Edit/duplicate gas cell", nextMenu=menuChooseLayerToEdit)
    plotLayerEntry = Entry("Plot gas cell", nextMenu=menuChoosePlotType)
    menuGasCell = Menu('Gas cell simulator', [createLayerEntry, editLayerEntry, plotLayerEntry], previousMenu=menuMain)
    return menuGasCell


def chooseDirection(profileList):
    lookUpEntry = Entry('Looking up', nextFunction=inputHeight, functionParams={'profiles': profileList,
                                                                                'direction': 'up'})
    lookDownEntry = Entry('Looking down', nextFunction=inputHeight, functionParams={'profiles': profileList,
                                                                                    'direction': 'down'})
    menuChooseDirection = Menu('Choose direction to look', [lookUpEntry, lookDownEntry], previousMenu=plotAtmTransferMenu)
    return menuChooseDirection


def chooseDirectionComponents(profile):
    lookUpEntry = Entry('Looking up', nextFunction=inputHeightComponents, functionParams={'profiles': profile,
                                                                                'direction': 'up'})
    lookDownEntry = Entry('Looking down', nextFunction=inputHeightComponents, functionParams={'profiles': profile,
                                                                                    'direction': 'down'})
    menuChooseDirection = Menu('Choose direction to look', [lookUpEntry, lookDownEntry],
                               previousMenu=plotProfileComponentsMenu)
    return menuChooseDirection


def plotProfileComponentsMenu(param=None):
    entryList = []
    profileList = util.getCompletedTransmissionList()
    for profile in profileList:
        entryList.append(Entry('%s' % profile, nextMenu=chooseDirectionComponents, functionParams=profile))
    menuAtmTransfer = Menu('Choose atmosphere', entryList, previousMenu=chooseAtmTransferBuildProfile, multiChoice=True)
    return menuAtmTransfer


def buildProfile(profileList):
    saveAbsData = pyradClasses.yesOrNo("Would you like to save abs coef data?\n"
                                       "This takes up quite a bit of space and generally isn't needed. %s"
                                       % util.limeText('(y/n)'))
    # initiate planet variable, in hopes of stopping memory leak
    planet = None

    for profile in profileList:
        print('Building %s on setting %s' % (profile.name, pyradClasses.settings.setting))
        planet = pyradClasses.createCustomPlanet(profile.name)
        overwrite = True
        progress, time = util.profileProgress(planet.folderPath)
        if util.profileComplete(planet.folderPath):
            overwrite = pyradClasses.yesOrNo("Data for this profile and setting seems to exist.\n"
                                             "Do you wish to overwrite it? %s" % util.limeText('(y/n)'))
        elif progress:
            resume = pyradClasses.yesOrNo('Partial data exists for this setting, do you wish to resume? %s' % util.limeText('(y/n)'))
            if not resume:
                overwrite = pyradClasses.yesOrNo("Are you sure, choosing yes will erase all previous data.\n"
                                                 "Do you wish to overwrite it? %s" % util.limeText('(y/n)'))
            else:
                planet.processLayers(verify=False, moleculeSpecific=True)
                overwrite = False
        if overwrite:
            util.emptyProfileDirectory(planet.folderPath)
            planet.processLayers(verify=False, moleculeSpecific=True)
        print('Creating transmission...')
        pyradClasses.processTransmissionBySingleLayer(planet.folderPath)
        if not saveAbsData:
            util.clearAbsData(planet.folderPath)
        planet.clearData()
    return chooseAtmTransferBuildProfile()


def buildProfilesMenu(param=None):
    entryList = []
    profileList = util.getPyrFileList()
    for profile in profileList:
        entryList.append(Entry('%s' % profile[:-4], nextFunction=buildProfile, functionParams=profile))
    menuBuildProfile = Menu('Which profile(s) do you want to build', entryList, previousMenu=chooseAtmTransferBuildProfile, multiChoice=True)
    return menuBuildProfile


def chooseAtmTransferBuildProfile(param=None):
    plotProfilesEntry = Entry('Plot profile(s)', nextMenu=plotAtmTransferMenu)
    plotProfileAndComponents = Entry('Plot single atm and components', nextMenu=plotProfileComponentsMenu)
    buildProfilesEntry = Entry('Build profile(s)', nextMenu=buildProfilesMenu)
    menuChooseProfileBuild = Menu('Plot or build profiles', [plotProfilesEntry, plotProfileAndComponents, buildProfilesEntry], previousMenu=menuMain)
    return menuChooseProfileBuild


def plotAtmTransferMenu(param=None):
    profileList = util.getCompletedTransmissionList()
    entryList = []
    for profile in profileList:
        entryList.append(Entry('%s' % profile, nextMenu=chooseDirection, functionParams=profile))
    menuAtmTransfer = Menu('Choose atmosphere', entryList, previousMenu=menuMain, multiChoice=True)
    return menuAtmTransfer


def chooseThemeMenu(originalMenu):
    themeList = pyradClasses.theme.listOfThemes
    entryList = []
    for theme in themeList:
        entryList.append(Entry('%s' % theme, nextFunction=loadTheme, functionParams=(theme, originalMenu)))
    menuChooseTheme = Menu('Choose theme', entryList, previousMenu=settingsMenu)
    return menuChooseTheme


def menuMain():
    gasCellEntry = Entry('Gas cell simulator', nextMenu=gasCellMenu)
    atmosphereTransferEntry = Entry('Atmosphere transmission', nextMenu=chooseAtmTransferBuildProfile)
    planckPlotEntry = (Entry('Plot planck curves', nextMenu=menuPlanckType))
    mainMenu = Menu('Main menu', [gasCellEntry, atmosphereTransferEntry, planckPlotEntry])
    return mainMenu


def settingsMenu(previousMenu=menuMain):
    lowSetting = Entry('low (intensity > 1e-21)', nextFunction=changeSettings, functionParams=('low', previousMenu))
    midSetting = Entry('mid (intensity > 1e-24)', nextFunction=changeSettings, functionParams=('mid', previousMenu))
    hiSetting = Entry('high (all absorption lines)', nextFunction=changeSettings, functionParams=('high', previousMenu))
    loadTheme = Entry('Change theme', nextMenu=chooseThemeMenu, functionParams=previousMenu)
    menuSettings = Menu('Choose level of detail', [lowSetting, midSetting, hiSetting, loadTheme], previousMenu=previousMenu)
    return menuSettings


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

